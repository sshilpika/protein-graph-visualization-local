# Visualizing Dynamically Generated Protein Interactions 

A Python library for visualizing dynamically generated protein interactions using 3D node-link layout.

[DEMO](./demo/Protein_Interaction_Demo.mp4)

## Installation

1. Requires `python=3.8`

   ```
   conda create -n myenv python=3.8
   conda activate myenv
   ```
2. Install packages
    ```
    cd src/
    pip install -r requirements.txt
    ```
Note: You may need to install [pygraphvis](https://pygraphviz.github.io/documentation/stable/install.html) using conda forge:
`conda install --channel conda-forge pygraphviz`

3. Start server
   ```
   ./start.sh
   ```

4. Open browser 
    ```angular2html
    http://127.0.0.1:5000/index
    ```