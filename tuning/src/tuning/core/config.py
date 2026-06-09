"""Load an experiment config from YAML.

A config just names the plugins to use, e.g.:

    component: clm
    sampler: history_matching
    runner: derecho

Look the classes up with get_component/get_sampler/get_runner (registry.py).
"""

import yaml


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)
