# Software_Example

An example that illustrates "good" and "bad" software design principles for a group meeting.
The examples use [`model.py`](model.py), a 2D QG model hosted [here](https://github.com/joernc/QGModel).

See [`example.ipynb`](example.ipynb) for a "bad" example and [`run.py`](run.py) and [`postprocess.py`](postprocess.py) for "good" examples, making use of the two "utils" modules.

Calkit can be installed with `uv tool install calkit-python`.

Then, to run the pipeline, which includes both the "good" and "bade" examples,
execute:

```sh
calkit run
```

Subsequent calls will skip stages if they or their inputs haven't changed.
It's also possible to open up the notebook and run it from top to bottom
using the uv-venv environment that is created in `.venv2`.
