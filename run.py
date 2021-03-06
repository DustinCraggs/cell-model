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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def run(args):
    with open(args.configuration) as f:
        config = json.load(f)
    _check_for_required_vars(config)

    overrides = []
    if 'variables' in config['experiment']:
        overrides = _get_overrides(config['experiment']['variables'])

    config_paths = [_get_config_path(config, o) for o in overrides]

    if args.overwrite:
        for path in config_paths:
            if os.path.exists(path):
                shutil.rmtree(os.path.dirname(path))

    for o in overrides:
        _write_configuration(copy.deepcopy(config), o)

    for path in config_paths:
        if args.synchronous:
            subprocess.run(['sbatch', '--wait', config['experiment']['job_script'], path])
            job_logs = glob.glob('slurm-*.out')
            print('Job complete; pausing for output')
            # Temporary: Wait for system to finish writing to log:
            time.sleep(6)
            subprocess.run(['cat'] + job_logs)
            subprocess.run(['mv'] + job_logs + [os.path.dirname(path)])
        else:
            subprocess.run(['sbatch', config['experiment']['job_script'], path])
            time.sleep(2)

    if not args.synchronous:
        time.sleep(6)
        job_logs = glob.glob('slurm-*.out')
        subprocess.run(['cat'] + job_logs)
        subprocess.run(['mv'] + job_logs + [os.path.dirname(os.path.dirname(path))])

def vis(args):
    with open(args.configuration) as f:
        config = json.load(f)
    _check_for_required_vars(config)

    overrides = []
    if 'variables' in config['experiment']:
        overrides = _get_overrides(config['experiment']['variables'])

    if 'plots' not in config['experiment']:
        sys.exit('No \'plots\' found in configuration')

    plots = config['experiment']['plots']
    if 'individual' in plots and plots['individual'] == True:
        for o in overrides:
            print('Plotting individual: {}'.format(o))
            path = _get_config_path(config, o)
            _plot_individual(path, overrides, config)
    print('Plotting runtime')
    if 'runtime' in plots and plots['runtime'] == True:
        _plot_runtime(config, overrides, args)

    print('Plotting metrics')
    if 'metrics' in plots and plots['metrics'] == True:
        _plot_metrics(config, overrides, args)

def combine(args):
    pass

def _plot_runtime(config, overrides, args):
    paths = [_get_config_path(config, o) for o in overrides]
    try:
        indep_var = config['experiment']['independent_variable']
    except KeyError:
        sys.exit('No \'independent_variable\' specified in experiment config')
    indep_vals = _get_indep_vals(indep_var, overrides)

    dfs = []
    for p in paths:
        directory = os.path.dirname(p)
        dfs.append(pd.read_csv(os.path.join(directory, config['output']['runtime']['file'])))
    times = pd.concat(dfs, axis=0)
    print(times)
    print(indep_vals)
    times['indep_var'] = indep_vals
    times['total_time'] = times['total_time']/1000
    times.plot(x='indep_var')
    basedir = os.path.dirname(os.path.dirname(paths[0]))
    plt.title('Simulation run time (1000 iterations)')
    plt.ylabel('Run time (s)')
    plt.xlabel(config['experiment']['independent_variable'])
    plt.savefig(os.path.join(basedir, 'runtime.png'), dpi=300)
    plt.close()

def _plot_metrics(config, overrides, args):
    print('plotting metrics')
    paths = [_get_config_path(config, o) for o in overrides]
    try:
        indep_var = config['experiment']['independent_variable']
    except KeyError:
        sys.exit('No \'independent_variable\' specified in experiment config')
    indep_vals = _get_indep_vals(indep_var, overrides)

    dfs = []
    for p in paths:
        directory = os.path.dirname(p)
        df = pd.read_csv(os.path.join(directory, config['output']['statistics']['file']))
        df = df.replace([np.inf, -np.inf], 0)
        dfs.append(pd.DataFrame(df.mean()).T)

    metrics = pd.concat(dfs, axis=0)
    basedir = os.path.dirname(os.path.dirname(paths[0]))
    metrics['indep_var'] = indep_vals
    metrics.to_csv(os.path.join(basedir, 'metrics.csv'))

    metrics.iloc[:, [1,-1]].plot(x='indep_var')
    plt.title('Average number of living cells')
    plt.ylabel('Average number of living cells')
    plt.xlabel(config['experiment']['independent_variable'])
    plt.savefig(os.path.join(basedir, 'number_of_cells.png'), dpi=300)
    plt.close()

    metrics.iloc[:, [2,3,4,-1]].plot(x='indep_var')
    plt.title('Cell resource summary')
    plt.ylabel('Average quantity of each resource')
    plt.xlabel(config['experiment']['independent_variable'])
    plt.savefig(os.path.join(basedir, 'cell_metrics.png'), dpi=300)
    plt.close()

    metrics.iloc[:, [5,6,-1]].plot(x='indep_var')
    plt.title('Environmental resource summary')
    plt.ylabel('Total quantity of each resource')
    plt.xlabel(config['experiment']['independent_variable'])
    plt.savefig(os.path.join(basedir, 'environment_metrics.png'), dpi=300)
    plt.close()

    metrics.iloc[:, [5,-1]].plot(x='indep_var')
    plt.title('Environmental chemicals')
    plt.ylabel('Total quantity of chemicals in environment')
    plt.xlabel(config['experiment']['independent_variable'])
    plt.savefig(os.path.join(basedir, 'environment_chem.png'), dpi=300)
    plt.close()

    metrics.iloc[:, [6,-1]].plot(x='indep_var')
    plt.title('Environmental waste')
    plt.ylabel('Total quantity of waste in environment')
    plt.xlabel(config['experiment']['independent_variable'])
    plt.savefig(os.path.join(basedir, 'environment_waste.png'), dpi=300)
    plt.close()

def _get_indep_vals(indep_var, overrides):
    vals = []
    for o in overrides:
        for var, val in o:
            if var == indep_var:
                vals.append(val)
    return vals

def _plot_individual(path, overrides, config):
    directory = os.path.dirname(path)
    stats = pd.read_csv(os.path.join(directory, config['output']['statistics']['file']))

    stats.iloc[:, 0:5].plot(x='iteration')
    plt.xlabel('Iteration number')
    plt.savefig(os.path.join(directory, 'CellStatistics.png'), dpi=300)
    plt.close()

    for col in stats.columns[1:]:
        stats[['iteration', col]].plot(x='iteration')
        plt.xlabel('Iteration number')
        plt.savefig(os.path.join(directory, '{}.png'.format(col)), dpi=300)
        plt.close()

    stats.iloc[:, [0, 5, 6]].plot(x='iteration')
    plt.xlabel('Iteration number')
    plt.savefig(os.path.join(directory, 'EnvironmentStatistics.png'), dpi=300)
    plt.close()

def _get_config_path(config, override):
    directory = config['experiment']['output_directory']
    filename = '_'.join(''.join(str(v) for v in t) for t in override)
    return os.path.join(directory, filename, 'config.json')

def _write_configuration(config, override):
    path = _get_config_path(config, override)
    output_dir = os.path.dirname(path)
    os.makedirs(output_dir, exist_ok=True)

    for var, val in override:
        config_var = config
        properties = [int(v) if v.isdigit() else v for v in var.split('.')]
        for v in properties[:-1]:
            config_var = config_var[v]
        config_var[properties[-1]] = val
    _prepend_output_directory(config, 'output.video.energy', output_dir)
    _prepend_output_directory(config, 'output.video.chemical', output_dir)
    _prepend_output_directory(config, 'output.video.toxin', output_dir)
    _prepend_output_directory(config, 'output.statistics.file', output_dir)
    _prepend_output_directory(config, 'output.runtime.file', output_dir)

    with open(path, 'w') as f:
        json.dump(config, f)
    return path

def _prepend_output_directory(config, prop_path, output_dir):
    props = prop_path.split('.')
    path_prop = config
    try:
        for prop in props[:-1]:
            path_prop = path_prop[prop]
        if len(path_prop[props[-1]]) > 0:
            path_prop[props[-1]] = os.path.join(output_dir, path_prop[props[-1]])
        # TODO: Check if path valid
    except KeyError:
        print('no : {}'.format(prop_path))
        return

def _get_overrides(variables):
    overrides = []
    for var, vals in variables.items():
        if var == 'iterations':
            overrides.extend(_merge_configs(o) for o in zip(*_get_overrides(vals)))
        elif var == 'permutations':
            overrides.extend(_merge_configs(o) for o in itertools.product(*_get_overrides(vals)))
        else:
            if not isinstance(vals, list):
                vals = [vals]
            overrides.append(list([tuple((var, v))] for v in vals))
    return overrides

def _merge_configs(configs):
    return list(itertools.chain.from_iterable(configs))

def _check_for_required_vars(config):
    if 'experiment' not in config:
        sys.exit('No \'experiment\' found in configuration')
    if 'output_directory' not in config['experiment']:
        sys.exit('No \'output_directory\' found in experiment configuration')
    if 'job_script' not in config['experiment']:
        sys.exit('No \'job_script\' found in experiment configuration')

def _parse_args():
    parser = argparse.ArgumentParser(description='Cell model run script and visualiser')
    subparsers = parser.add_subparsers(title="Subcommands", dest='subcommand', required=True)

    # Run parser:
    run_parser = subparsers.add_parser('run', help='Run an experiment')
    run_parser.add_argument('configuration', help='Path to a JSON configuration file')
    run_parser.add_argument('-s', '--synchronous', action='store_true', help='Run one phoenix job at a time')
    run_parser.add_argument('-o', '--overwrite', action='store_true',
        help='Overwrite previous results for this experiment')
    run_parser.set_defaults(func=run)
    
    # Visualisation parser:
    vis_parser = subparsers.add_parser('vis', help='Visualise results of an experiment')
    vis_parser.add_argument('configuration', help='Path to a JSON configuration file')
    vis_parser.set_defaults(func=vis)

    return parser.parse_args()

if __name__=='__main__':
    args = _parse_args()
    args.func(args)