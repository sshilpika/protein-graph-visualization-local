from visg import data_path, master_file, data_part_width, master_filename

import re
import os
from os import listdir
from os.path import isfile, join
import json
import subprocess
import pygraphviz
from pygraphviz import AGraph
import networkx as nx
from networkx.readwrite import json_graph
from joblib import Parallel, delayed
from operator import itemgetter

class Protein_Graph:
    """ Class for reading and updating protein to protein interactions in a graph format.
        The dot file is read and data partitions (specified in data_part_width)
        are stored in json files to be read by the frontend.
        finalK : count of files to be visualized
        processed_line_counts : number of lines already processed from the master file
        minlink_count, maxlink_count : minimum and maximum number of links to be visualized
    """

    data_part_width = data_part_width
    data_path = data_path
    clean_master_file = None
    finalK = 0
    minlink_count = 0
    maxlink_count = 3000000
    processed_line_counts = 0


    def __init__(self):
        self.A = None
        self.G = None

    def clean_dot_file(self, master_file_name):

        f = open(os.path.join(Protein_Graph.data_path, master_file_name))
        s=f.read().replace("-", "_" ).replace("_>", "->")
        f.close()

        self.clean_master_file = "clean_"+master_file_name

        with open(os.path.join(Protein_Graph.data_path, self.clean_master_file), 'w') as f:
            f.write(s)

        return self.clean_master_file


    def read_graph_from_file(self, master_file_name):
        if "temp" in master_file_name:
            Protein_Graph.processed_line_counts += Protein_Graph.wc_l(Protein_Graph.data_path, master_file_name)
            Protein_Graph.processed_line_counts -= 2
        else:
            Protein_Graph.processed_line_counts = Protein_Graph.wc_l(Protein_Graph.data_path, master_file_name)
        self.A = AGraph(string=open(os.path.join(Protein_Graph.data_path,master_file_name)).read())
        self.G = nx.DiGraph(self.A)


    def remove_nodes_from_graph(self):
        remove = [node for node,degree in dict(self.G.degree()).items() if degree < Protein_Graph.minlink_count or degree > Protein_Graph.maxlink_count]
        self.G.remove_nodes_from(remove)
        remove = [node for node,degree in dict(self.G.degree()).items() if degree == 0 ]
        self.G.remove_nodes_from(remove)

    @staticmethod
    def wc_l(directory_path, filename):
        line_count = int(subprocess.check_output(["wc", "-l", os.path.join(directory_path, filename)]).split()[0])
        return line_count

    # format graph and save in json files
    @staticmethod
    def get_graph(minlink_count = 0, maxlink_count = 3000000, master_file_name= master_filename, remove_nodes = True, set_iteratons = True):

        pgraph = Protein_Graph()
        Protein_Graph.minlink_count = minlink_count
        Protein_Graph.maxlink_count = maxlink_count
        clean_master_file_name = pgraph.clean_dot_file(master_file_name)
        pgraph.read_graph_from_file(clean_master_file_name)

        if remove_nodes:
            pgraph.remove_nodes_from_graph()


        dict_json_ = json_graph.node_link_data(pgraph.G)
        data_main = {}
        data_main["nodes"] = dict_json_["nodes"]
        data_main["links"] = dict_json_["links"]
        for i, l in enumerate(data_main['links']):
            l.update({'id': i})
        print("Data Cleaning Done...")

        step = Protein_Graph.data_part_width

        if set_iteratons:
            Protein_Graph.finalK = 0

        # joblib processing with threads
        def process_graph(k, i, step, data_main, G, data_path, data_part_width, finalK):

            data_main_part = {}
            data_main_part["nodes"] = data_main["nodes"][ i-step : i ]
            all_nodes = [d["id"] for d in data_main_part["nodes"]]

            for n in data_main_part["nodes"]:
                n["neighbors"] = [nei for nei in list(nx.all_neighbors(G, n["id"])) if nei in all_nodes]
                n["links"] = []
                n['node_neighbor_count'] = 0

                for s, t in G.in_edges(n["id"]):
                    if s in all_nodes and t in all_nodes:
                        n_in = [0] if (s == n["id"]) else [s1 for s1 in list(nx.all_neighbors(G, s)) if s1 != n["id"] and s1 in all_nodes ]
                        n_out = [0] if (t == n["id"]) else [t1 for t1 in list(nx.all_neighbors(G, t)) if t1 != n["id"] and t1 in all_nodes]
                        n['node_neighbor_count'] = n['node_neighbor_count'] + len(set(n_in+n_out))
                        link = {
                            "source":{"id":s}, "target":{"id":t},
                            "source_neighbors": len(n_in),
                            "target_neighbors": len(n_out)
                        }
                        link.update(G.get_edge_data(s, t))
                        n["links"].append(link)

                for s, t in G.out_edges(n["id"]):
                    if s in all_nodes and t in all_nodes:
                        n_in = [0] if (s == n["id"]) else [s1 for s1 in list(nx.all_neighbors(G, s)) if s1 != n["id"] and s1 in all_nodes ]
                        n_out = [0] if (t == n["id"]) else [t1 for t1 in list(nx.all_neighbors(G, t)) if t1 != n["id"] and t1 in all_nodes ]
                        n['node_neighbor_count'] = n['node_neighbor_count'] + len(set(n_in+n_out))
                        link = {
                            "source":{"id":s}, "target":{"id":t},
                            "source_neighbors": len(n_in),
                            "target_neighbors": len(n_out)
                        }
                        link.update(G.get_edge_data(s, t))
                        n["links"].append(link)


            data_main_part["links"] = [link for link in data_main["links"] if link["source"] in all_nodes and link["target"] in all_nodes]


            with open(os.path.join(data_path,"graph_master_part"+str(data_part_width)+"_"+str(k + finalK)+".json"), 'w') as f:
                json.dump(data_main_part, f, ensure_ascii=False)

            return (step, k)

        results = Parallel(n_jobs=10)(delayed(process_graph)(k, i, step, data_main, pgraph.G, Protein_Graph.data_path, Protein_Graph.data_part_width, Protein_Graph.finalK) for k,i in enumerate(range(step,len(data_main["nodes"])+step, step)))
        increment_k  = max(results, key=itemgetter(1))[1] + 1
        print(results, increment_k)

        Protein_Graph.finalK = increment_k if set_iteratons else increment_k + Protein_Graph.finalK
        return Protein_Graph.finalK