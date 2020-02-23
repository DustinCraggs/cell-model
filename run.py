import argparse
import subprocess
import json
import itertools

import pandas as pd
import matplotlib.pyplot as plt

def run(args):
    with open(args.configuration) as f:
        config = json.load(f)
    _check_for_required_vars(config)

    overrides = _get_overrides(config['experiment']['variables'])


def vis(args):
    print(json.dumps(config, indent=4))

def combine(args):
    pass

def _get_overrides(variables):
    overrides = []
    for var, vals in variables.items():
        if var == 'iterations':
            overrides.extend(list(o for o in zip(*_get_overrides(vals))))
        elif var == 'permutations':
            overrides.extend(list(list(o) for o in itertools.product(*_get_overrides(vals))))
            # itertools.permutations(_get_overrides(val))
        else:
            if not isinstance(vals, list):
                vals = [vals]
            overrides.append(list([tuple((var, v))] for v in vals))
    print(overrides)
    return overrides

def _check_for_required_vars(config):
    if 'experiment' not in config:
        sys.exit('No \'experiment\' found in configuration')
    if 'variables' not in config['experiment']:
        sys.exit('No \'variables\' found in experiment configuration')


def _parse_args():
    parser = argparse.ArgumentParser(description='Cell model run script and visualiser')
    subparsers = parser.add_subparsers(title="Subcommands", dest='subcommand', required=True)

    # Run parser:
    run_parser = subparsers.add_parser('run', help='Run an experiment')
    run_parser.add_argument('configuration', help='Path to a JSON configuration file')
    run_parser.set_defaults(func=run)
    
    # Visualisation parser:
    vis_parser = subparsers.add_parser('vis', help='Visualise results of an experiment')
    vis_parser.add_argument('configuration', help='Path to a JSON configuration file')
    vis_parser.set_defaults(func=run)

    return parser.parse_args()

if __name__=='__main__':
    args = _parse_args()
    args.func(args)