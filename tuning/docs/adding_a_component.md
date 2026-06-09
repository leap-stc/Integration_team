# Adding a component

A component is a CESM model. It lists its tunable parameters, writes them into a
case, and turns a run into diagnostics. Subclass `Component`, implement three
methods, and register a name.

```python
from tuning.core.interfaces import Component
from tuning.core.parameters import Parameter, ParameterSet
from tuning.core.registry import register_component


@register_component("my_model")
class MyModel(Component):
    def parameters(self):
        return ParameterSet([Parameter("x", low=0.0, high=10.0, default=5.0)])

    def apply(self, case_dir, params):
        # write params into the CESM case (namelist / xml / paramfile)
        ...

    def compute_metrics(self, run_output):
        # return a numpy array matching your Observations.targets
        ...
```

Put it in `src/tuning/components/my_model.py` and import it in
`components/__init__.py`. See `components/clm.py` for an example.
