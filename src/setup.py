# from setuptools import setup, find_packages
from distutils.core import setup
setup(
    name="visg",
    version = "0.1",
    packages=["flask", "Werkzeug", "Flask-Cors", "Flask-Sijax", "pygraphviz", "networkx", "joblib", "watchdog"],
)