# mpr-python
Jupyter-notebook for Muti-scale Parameter Regionalization




### Directory structure

1 ./notebook

notebook to explore mpr

2 ./mpr

a collection of custom modules


### To use mpr modules:

Install the mpr by

```bash
cd mpr-python 
(optional) conda activate $ENVIRONMENT_NAME
pip install -e .
```
after that, mpr module can be imported

```python
import mpr 
```

or

if you would like not to install mpr module,
```python
insert sys.path.append(<PATH TO mpr-python>)
```
e.g., calling from ./notebook/xx.ipyb, insert
```python
sys.path.append('../')
```
