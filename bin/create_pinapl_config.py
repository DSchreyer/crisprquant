#!/usr/bin/env python
import argparse
import yaml


def read_param() -> dict:
    """
    read parameters from terminal
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ScreenType", help="type of screen ['enrichment'/'depletion']", type=str, choices=["enrichment", "depletion"]
    )
    parser.add_argument("--LibFilename", help="filename of library spreadsheet", type=str)
    parser.add_argument("--seq_5_end", help="5'-sequence adapter", type=str)
    parser.add_argument("--CtrlPrefix", help="Name of control", type=str)
    parser.add_argument("--NonTargetPrefix", help="prefix for non-targeting sgRNAs in library", type=str)
    parser.add_argument("--NumGuidesPerGene", help="number of sgRNAs per gene", type=int)

    args = parser.parse_args()
    # create a dictionary to store arguments
    args_dict = dict()
    for arg in vars(args):
        args_dict[arg] = getattr(args, arg)

    return args_dict


def yaml_converter(parameters: dict) -> yaml:
    """
    converts given parameters to a yaml file
    """
    # TODO: rename testparse.yaml to configuration.yaml
    with open("testparse.yaml", "w") as file:
        documents = yaml.dump(parameters, file)


def main():
    """
    main function
    """
    # create config file
    yaml_converter(parameters=read_param())


if __name__ == "__main__":
    main()
