# Software_Example

An example that illustrates "good" and "bad" software design principles for a group meeting.
The examples use [`model.py`](model.py), a 2D QG model hosted [here](https://github.com/joernc/QGModel).

See [`example.ipynb`](example.ipynb) for a "bad" example and [`run.py`](run.py) and [`postprocess.py`](postprocess.py) for "good" examples, making use of the two "utils" modules.

Calkit can be installed with `uv tool install calkit-python`.

Then, to run the pipeline, which includes both the "good" and "bad" examples,
execute:

```sh
calkit run
```

Subsequent calls will skip stages if they or their inputs haven't changed.
It's also possible to open up the notebook in JupyterLab
using its dedicated environment (created in `.venv2`) with:

```sh
calkit xenv --name nb jupyter lab
```

If the pipeline ran successfully, you'll notice the cell that runs the
simulation will exit quickly since its results have been cached.
It's possible to push this cache to the cloud so collaborators will also
not need to rerun the simulation cell.
