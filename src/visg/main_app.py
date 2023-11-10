from visg import app
from flask import render_template
from flask import Flask, g
import flask_sijax
from visg import data_path, master_file, data_part_width, master_filename, new_data_master_filename
# from visg.scripts.listener import Listener
from visg.scripts.graph_processor import Protein_Graph

import re
import os
from os import listdir
from os.path import isfile, join
import json
import pygraphviz
from pygraphviz import AGraph
import networkx as nx
from networkx.readwrite import json_graph
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import os
import itertools



class Listener:
    """ Class for watching master file updates.
        When the master file has been updated, the new data is stored in a temp dot file
        and the number of processed lines is updated.
    """
    DIRECTORY_TO_WATCH = data_path
    FILE_TO_WATCH = master_filename
    FILE_TO_WRITE = new_data_master_filename
    watchFlag = True

    def __init__(self):
        self.observer = Observer()

    def run(self):
        event_handler = Handler()
        self.observer.schedule(event_handler, self.DIRECTORY_TO_WATCH, recursive=True)
        self.observer.start()
        try:
            while True:
                time.sleep(5)
        except:
            self.observer.stop()
        self.observer.join()


class Handler(FileSystemEventHandler):
    @staticmethod
    def on_modified(event):
        if event.is_directory:
            return None
        elif event.src_path.endswith(Listener.FILE_TO_WATCH) and Listener.watchFlag:
            lineC = Protein_Graph.wc_l(Listener.DIRECTORY_TO_WATCH, Listener.FILE_TO_WATCH)
            processed_line_counts = Protein_Graph.processed_line_counts
            print("file changed", processed_line_counts, Listener.FILE_TO_WATCH, Listener.watchFlag, lineC)
            if processed_line_counts < lineC:
                lines = []
                with open(os.path.join(Listener.DIRECTORY_TO_WATCH, Listener.FILE_TO_WATCH), "r") as f:
                    a = processed_line_counts-1
                    a = (abs(a)+a)//2
                    lines = [line for line in itertools.islice(f, a, lineC)]

                with open(os.path.join(Listener.DIRECTORY_TO_WATCH, Listener.FILE_TO_WRITE), "w") as f:
                    f.write("digraph G {\n")
                    for line in lines:
                        f.write(line)

                    processed_line_counts = lineC
                finalK = Protein_Graph.get_graph(Protein_Graph.minlink_count, Protein_Graph.maxlink_count, Listener.FILE_TO_WRITE, True, False)
                print(f'new final counter value is {finalK}')



def toggle_listener(watch=None):
    if watch == None:
        Listener.watchFlag = not Listener.watchFlag
    else:
        Listener.watchFlag = watch

def read_from_dot(dot_file):
    A = AGraph(string=open(dot_file).read())
    G = nx.DiGraph(A)
    dict_json_ = json_graph.node_link_data(G)

def get_graph_partions(obj_response, filename, reset):
    f = open(os.path.join(data_path,"graph_master_part"+str(data_part_width)+"_"+str(filename)+".json"))
    data_partition = json.load(f)
    obj_response.script("addGraphData("+str(data_partition)+","+reset+")")

def set_nodelink_limit(obj_response, minlink_count, maxlink_count):
    toggle_listener(False)
    finalK = Protein_Graph.get_graph(minlink_count, maxlink_count, master_filename, True, True)
    toggle_listener(True)
    obj_response.script("reloadGraphData(reset = true, stopAt = "+str(finalK)+")")

def startListener(obj_response):
    my_listener = Listener()
    my_listener.run()

def setWatchFlag(obj_response):
    print(f'BBBWatch file flag = {Listener.watchFlag}')
    toggle_listener()
    print(f'WWWatch file flag = {Listener.watchFlag}')

def checkFinalK(obj_response, current_stop_at):
    obj_response.script("updateStopAt("+str(Protein_Graph.finalK)+")")


@app.route('/')
def hello():
    return render_template("index_main.html")

@flask_sijax.route(app, '/index')
def index():

    if g.sijax.is_sijax_request:
        g.sijax.register_callback('getDataPartions', get_graph_partions)
        g.sijax.register_callback('setNodeLinkLimit', set_nodelink_limit)
        g.sijax.register_callback('watchMasterFile', startListener)
        g.sijax.register_callback('setWatchFlag', setWatchFlag)
        g.sijax.register_callback('checkStopAt', checkFinalK)
        return g.sijax.process_request()

    return render_template("index_main_3D.html")

