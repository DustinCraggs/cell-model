import argparse
import subprocess
import json
import itertools
import os
import sys
import glob
import shutil
import copy
import time

from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

maxGenomeNum = 10
numberOfVariablesBeforeGenomes = 8

def run(args):
    with open(args.configuration) as f:
        config = json.load(f)
    _check_for_required_vars(config)

    overrides = []
    if "variables" in config["experiment"]:
        overrides = _get_overrides(config["experiment"]["variables"])

    # Only supports one experiment:
    overrides = overrides[0]

    config_paths = [_get_config_path(config, o) for o in overrides]

    if args.overwrite:
        for path in config_paths:
            if os.path.exists(path):
                shutil.rmtree(os.path.dirname(path))

    for o in overrides:
        _write_configuration(copy.deepcopy(config), o)


    print(config_paths)

    for path in config_paths:
        if not args.phoenix:
            subprocess.run(["./build/simulate", path])
        elif args.synchronous:
            subprocess.run(
                ["sbatch", "--wait", config["experiment"]["job_script"], path]
            )
            job_logs = glob.glob("slurm-*.out")
            print("Job complete; pausing for output")
            # Temporary: Wait for system to finish writing to log:
            time.sleep(6)
            subprocess.run(["cat"] + job_logs)
            subprocess.run(["mv"] + job_logs + [os.path.dirname(path)])
        else:
            subprocess.run(["sbatch", config["experiment"]["job_script"], path])
            time.sleep(1)


def vis(args):
    with open(args.configuration) as f:
        config = json.load(f)
    _check_for_required_vars(config)

    overrides = []
    if "variables" in config["experiment"]:
        overrides = _get_overrides(config["experiment"]["variables"])[0]

    if "plots" not in config["experiment"]:
        sys.exit("No 'plots' found in configuration")

    plots = config["experiment"]["plots"]
    genomeNum = config["model"]["genomeNum"];

    count = 1;
    currentGenomeNum = 1;
    length = len(overrides)
    seeds = length/genomeNum

    if "individual" in plots and plots["individual"] == True:
        for o in overrides:
            print("Plotting individual: {}".format(o))
            path = _get_config_path(config, o)
            _plot_individual(path, overrides, config, currentGenomeNum)
            if count%seeds == 0:
                currentGenomeNum+=1
            count+=1
    print("Plotting runtime")
    if "runtime" in plots and plots["runtime"] == True:
        _plot_runtime(config, overrides, args)

    print("Plotting metrics")
    if "metrics" in plots and plots["metrics"] == True:
        _plot_metrics(config, overrides, args)


def combine(args):
    pass


def _plot_runtime(config, overrides, args):
    paths = [_get_config_path(config, o) for o in overrides]
    try:
        indep_var = config["experiment"]["independent_variable"]
    except KeyError:
        sys.exit("No 'independent_variable' specified in experiment config")
    indep_vals = _get_indep_vals(indep_var, overrides)

    totalRuntime = 0

    dfs_by_indep_val = defaultdict(list)

    for x, p in zip(indep_vals, paths):
        directory = os.path.dirname(p)
        df = pd.read_csv(os.path.join(directory, config["output"]["runtime"]["file"]))
        dfs_by_indep_val[x].append(df)
        totalRuntime += df["total_time"]

    print(totalRuntime)

    basedir = os.path.dirname(os.path.dirname(paths[0]))
    totalRuntime.to_csv(os.path.join(basedir, "totalRuntime.csv"))

    dfs = []
    for indep_val, df_list in dfs_by_indep_val.items():
        df = pd.concat(df_list, axis=0)
        df["indep_var"] = indep_val
        dfs.append(pd.DataFrame(df.mean()).T)

    times = pd.concat(dfs, axis=0)

    print(times)
    print(indep_vals)

    times["total_time"] = times["total_time"] / 1000
    times.plot(x="indep_var")
    basedir = os.path.dirname(os.path.dirname(paths[0]))
    plt.title("Simulation run time (1000 iterations)")
    plt.ylabel("Run time (s)")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.savefig(os.path.join(basedir, "runtime.png"), dpi=300)
    plt.close()


def _plot_metrics(config, overrides, args):
    print("plotting metrics")
    paths = [_get_config_path(config, o) for o in overrides]
    print(paths);
    try:
        indep_var = config["experiment"]["independent_variable"]
    except KeyError:
        sys.exit("No 'independent_variable' specified in experiment config")
    indep_vals = _get_indep_vals(indep_var, overrides)

    dfs_by_indep_val = defaultdict(list)
    for x, p in zip(indep_vals, paths):
        directory = os.path.dirname(p)
        df = pd.read_csv(
            os.path.join(directory, config["output"]["statistics"]["file"])
        )
        df = df.replace([np.inf, -np.inf], 0)
        dfs_by_indep_val[x].append(df)

    dfs = []
    for indep_val, df_list in dfs_by_indep_val.items():
        df = pd.concat(df_list, axis=0)
        df["indep_var"] = indep_val
        dfs.append(pd.DataFrame(df.mean()).T)

    metrics = pd.concat(dfs, axis=0)
    basedir = os.path.dirname(os.path.dirname(paths[0]))
    os.makedirs(basedir, exist_ok=True)
    metrics.to_csv(os.path.join(basedir, "metrics.csv"))

    genomeNum = config["model"]["genomeNum"];

    # MAKE LEGEND MODULAR

    averageNumGenome = [];

    for x in range(genomeNum):
        averageNumGenome.append(numberOfVariablesBeforeGenomes+x)

    averageNumGenome.append(-1)

    metrics.iloc[:, averageNumGenome].plot(x="indep_var")
    plt.title("Average number of living cells by genome")
    plt.ylabel("Average living cells")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.legend(["genome1", "genome2", "genome3", "genome4", "genome5", "genome6", "genome7", "genome8", "genome9", "genome10"], loc = "upper left")
    plt.savefig(os.path.join(basedir, "genome_number_of_cells.png"), dpi=300)
    plt.close()

    averageEnergyGenome = []

    for x in range(genomeNum):
        averageEnergyGenome.append(numberOfVariablesBeforeGenomes+maxGenomeNum+x)

    averageEnergyGenome.append(-1)

    metrics.iloc[:, averageEnergyGenome].plot(x="indep_var")
    plt.title("Average number of cell energy by genome")
    plt.ylabel("Average cell energy")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.legend(["genome1", "genome2", "genome3", "genome4", "genome5", "genome6", "genome7", "genome8", "genome9", "genome10"], loc = "upper left")
    plt.savefig(os.path.join(basedir, "genome_average_energy.png"), dpi=300)
    plt.close()

    averageChemGenome = []

    for x in range(genomeNum):
        averageChemGenome.append(numberOfVariablesBeforeGenomes+(maxGenomeNum*2)+x)

    averageChemGenome.append(-1)

    metrics.iloc[:, averageChemGenome].plot(x="indep_var")
    plt.title("Average number of cell chemicals by genome")
    plt.ylabel("Average cell chemicals")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.legend(["genome1", "genome2", "genome3", "genome4", "genome5", "genome6", "genome7", "genome8", "genome9", "genome10"], loc = "upper left")
    plt.savefig(os.path.join(basedir, "genome_average_chemicals.png"), dpi=300)
    plt.close()

    metrics.iloc[:, [1, -1]].plot(x="indep_var")
    plt.title("Average number of living cells")
    plt.ylabel("Average number of living cells")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.savefig(os.path.join(basedir, "number_of_cells.png"), dpi=300)
    plt.close()

    metrics.iloc[:, [2, -1]].plot(x="indep_var")
    plt.title("Average cell size")
    plt.ylabel("Average cell size")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.savefig(os.path.join(basedir, "cell_size.png"), dpi=300)
    plt.close()

    metrics.iloc[:, [3, 4, 5, -1]].plot(x="indep_var")
    plt.title("Cell resource summary")
    plt.ylabel("Average quantity of each resource")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.savefig(os.path.join(basedir, "cell_metrics.png"), dpi=300)
    plt.close()

    metrics.iloc[:, [6, 7, -1]].plot(x="indep_var")
    plt.title("Environmental resource summary")
    plt.ylabel("Total quantity of each resource")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.savefig(os.path.join(basedir, "environment_metrics.png"), dpi=300)
    plt.close()

    metrics.iloc[:, [6, -1]].plot(x="indep_var")
    plt.title("Environmental chemicals")
    plt.ylabel("Total quantity of chemicals in environment")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.savefig(os.path.join(basedir, "environment_chem.png"), dpi=300)
    plt.close()

    metrics.iloc[:, [7, -1]].plot(x="indep_var")
    plt.title("Environmental waste")
    plt.ylabel("Total quantity of waste in environment")
    plt.xlabel(config["experiment"]["independent_variable"])
    plt.savefig(os.path.join(basedir, "environment_waste.png"), dpi=300)
    plt.close()


def _get_indep_vals(indep_var, overrides):
    vals = []
    for o in overrides:
        for var, val in o:
            if var == indep_var:
                vals.append(val)
    return vals


def _plot_individual(path, overrides, config, genomeNum):

    directory = os.path.dirname(path)
    stats = pd.read_csv(os.path.join(directory, config["output"]["statistics"]["file"]))

    stats.iloc[:, 0:5].plot(x="iteration")
    plt.xlabel("Iteration number")
    plt.savefig(os.path.join(directory, "CellStatistics.png"), dpi=300)
    plt.close()

    for col in stats.columns[1:]:
        stats[["iteration", col]].plot(x="iteration")
        plt.xlabel("Iteration number")
        plt.savefig(os.path.join(directory, "{}.png".format(col)), dpi=300)
        plt.close()

    stats.iloc[:, [0, 5, 6]].plot(x="iteration")
    plt.xlabel("Iteration number")
    plt.savefig(os.path.join(directory, "EnvironmentStatistics.png"), dpi=300)
    plt.close()

    cellNum = [0]
    legend = []

    for x in range(genomeNum):
        cellNum.append(numberOfVariablesBeforeGenomes+x)

    for x in range(genomeNum):
        legend.append("genome" + str(x+1))

    stats.iloc[:, cellNum].plot(x="iteration")
    plt.title("Number of cells by genome")
    plt.ylabel("Number of cells")
    plt.xlabel("iterations")
    plt.legend(legend, loc = "upper left")
    plt.savefig(os.path.join(directory, "genome_cell_num.png"), dpi=300)
    plt.close()

    cellEnergy = [0]

    for x in range(genomeNum):
        cellEnergy.append(numberOfVariablesBeforeGenomes+maxGenomeNum+x)

    stats.iloc[:, cellEnergy].plot(x="iteration")
    plt.title("Average cell energy by genome")
    plt.ylabel("Average energy")
    plt.xlabel("iterations")
    plt.legend(legend, loc = "upper left")
    plt.savefig(os.path.join(directory, "genome_average_energy.png"), dpi=300)
    plt.close()

    cellChem = [0]

    for x in range(genomeNum):
        cellChem.append(numberOfVariablesBeforeGenomes+maxGenomeNum*2+x)

    stats.iloc[:, cellChem].plot(x="iteration")
    plt.title("Average cell chemicals by genome")
    plt.ylabel("Average chemicals")
    plt.xlabel("iterations")
    plt.legend(legend, loc = "upper left")
    plt.savefig(os.path.join(directory, "genome_average_chem.png"), dpi=300)
    plt.close()

def _get_config_path(config, override):
    directory = config["experiment"]["output_directory"]
    filename = "_".join("".join(str(v) for v in t) for t in override)
    return os.path.join(directory, filename, "config.json")


def _write_configuration(config, override):
    path = _get_config_path(config, override)
    output_dir = os.path.dirname(path)
    os.makedirs(output_dir, exist_ok=True)

    for var, val in override:
        config_var = config
        properties = [int(v) if v.isdigit() else v for v in var.split(".")]
        for v in properties[:-1]:
            config_var = config_var[v]
        config_var[properties[-1]] = val
    _prepend_output_directory(config, "output.video.energy", output_dir)
    _prepend_output_directory(config, "output.video.chemical", output_dir)
    _prepend_output_directory(config, "output.video.toxin", output_dir)
    _prepend_output_directory(config, "output.video.genome", output_dir)
    _prepend_output_directory(config, "output.statistics.file", output_dir)
    _prepend_output_directory(config, "output.runtime.file", output_dir)

    with open(path, "w") as f:
        json.dump(config, f)
    return path


def _prepend_output_directory(config, prop_path, output_dir):
    props = prop_path.split(".")
    path_prop = config
    try:
        for prop in props[:-1]:
            path_prop = path_prop[prop]
        if len(path_prop[props[-1]]) > 0:
            path_prop[props[-1]] = os.path.join(output_dir, path_prop[props[-1]])
        # TODO: Check if path valid
    except KeyError:
        print("no : {}".format(prop_path))
        return


def _get_overrides(variables):
    overrides = []
    for var, vals in variables.items():
        if var.startswith("iterations"):
            overrides.append([_merge_configs(o) for o in zip(*_get_overrides(vals))])
        elif var.startswith("permutations"):
            overrides.append(
                [_merge_configs(o) for o in itertools.product(*_get_overrides(vals))]
            )
        else:
            if not isinstance(vals, list):
                vals = [vals]
            overrides.append(list([tuple((var, v))] for v in vals))
    return overrides


def _merge_configs(configs):
    return list(itertools.chain.from_iterable(configs))


def _check_for_required_vars(config):
    if "experiment" not in config:
        sys.exit("No 'experiment' found in configuration")
    if "output_directory" not in config["experiment"]:
        sys.exit("No 'output_directory' found in experiment configuration")
    if "job_script" not in config["experiment"]:
        sys.exit("No 'job_script' found in experiment configuration")


def _parse_args():
    parser = argparse.ArgumentParser(description="Cell model run script and visualiser")
    subparsers = parser.add_subparsers(
        title="Subcommands", dest="subcommand", required=True
    )

    # Run parser:
    run_parser = subparsers.add_parser("run", help="Run an experiment")
    run_parser.add_argument("configuration", help="Path to a JSON configuration file")
    run_parser.add_argument(
        "-p", "--phoenix", action="store_true", help="Run one phoenix job at a time"
    )
    run_parser.add_argument(
        "-s", "--synchronous", action="store_true", help="Run one phoenix job at a time"
    )
    run_parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="Overwrite previous results for this experiment",
    )
    run_parser.set_defaults(func=run)

    # Visualisation parser:
    vis_parser = subparsers.add_parser("vis", help="Visualise results of an experiment")
    vis_parser.add_argument("configuration", help="Path to a JSON configuration file")
    vis_parser.set_defaults(func=vis)

    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    args.func(args)
