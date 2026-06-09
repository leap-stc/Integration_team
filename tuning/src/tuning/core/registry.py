"""A name -> class lookup so configs can pick plugins by name.

Each plugin registers itself, e.g.:

    @register_sampler("random")
    class RandomSampler(Sampler): ...

Then `get_sampler("random")` returns the class.
"""

_samplers: dict[str, type] = {}
_components: dict[str, type] = {}
_runners: dict[str, type] = {}


def register_sampler(name):
    def add(cls):
        _samplers[name] = cls
        return cls
    return add


def register_component(name):
    def add(cls):
        _components[name] = cls
        return cls
    return add


def register_runner(name):
    def add(cls):
        _runners[name] = cls
        return cls
    return add


def get_sampler(name):
    return _samplers[name]


def get_component(name):
    return _components[name]


def get_runner(name):
    return _runners[name]
