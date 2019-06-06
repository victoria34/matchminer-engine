##
# imports
import sys
import os
import networkx as nx

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen

## constants.
FILE_URL = "http://oncotree.mskcc.org/api/tumor_types.txt"

## functions

def build_oncotree(file_path=False, metamaintype_col=7):

    # load the file.
    if not file_path:

        # fetch from inter-webs.
        #req = urllib2.Request(FILE_URL)
        response = urlopen(FILE_URL)
        the_page = response.read().decode('utf-8')

        # split into array.
        lines = the_page.strip().split("\n")

    else:

        # just open the file.
        with open(file_path, "r") as fin:
          lines = fin.readlines()

    # create a graph.
    g = nx.DiGraph()

    # add root node.
    g.add_node("root", text="root")
    root = "root"

    # parse the file.
    line_cnt = 0
    old_style = False
    for line in lines:

        # tokenize.
        tokens = line.strip().split("\t")

        # skip header.
        if line_cnt == 0:
            if len(tokens) == 9:
              old_style = True
            line_cnt += 1
            continue

        def try_to_set(i, tokens):
          try:
            v = tokens[i]
          except:
            v = None
          return v

        # check for version.
        history = try_to_set(11, tokens)
        metaumls = try_to_set(10, tokens)
        metanci = try_to_set(9, tokens)
        metacolor = try_to_set(8, tokens)
        metamaintype = try_to_set(metamaintype_col, tokens)
        level7 = try_to_set(6, tokens)
        level6 = try_to_set(5, tokens)
        level5 = try_to_set(4, tokens)
        level4 = try_to_set(3, tokens)
        level3 = try_to_set(2, tokens)
        level2 = try_to_set(1, tokens)
        level1 = try_to_set(0, tokens)

        if old_style:
          levels = [level1, level2, level3, level4, level5]
        else:
          levels = [level1, level2, level3, level4, level5, level6, level7]

        #print(metanci, metamaintype)
        #print(levels)
        #break

        # set root node.
        prev_n = root

        # build nodes all the way down.
        for i in range(0, len(levels)):

            # check if empty.
            if levels[i] == "" or levels[i] is None or len(levels[i].split('(')) < 2:
                continue

            # split into two.
            tmp = levels[i].split("(")
            val = tmp[0].strip().replace('"', '').replace("'", '')
            key = tmp[1].strip().replace("(", "").replace(")", "").replace('"', '').replace("'", '')

            # build node.
            g.add_node(key,
                       text=val,
                       metamaintype=metamaintype,
                       metacolor=metacolor,
                       metanci=metanci,
                       metaumls=metaumls,
                       history=history)
            n = key

            # add edge.
            g.add_edge(prev_n, n)

            # update previous node.
            prev_n = n

        # increment line count.
        line_cnt += 1

    # return the graph.
    return g


def get_basal(g, source):

    # get ancestors.
    nlist = nx.ancestors(g, source)

    # find the one with predecessor root.
    hit = None
    for n in nlist:

        # get predecessor.
        preds = list(g.predecessors(n))
        if len(preds) > 0 and preds[0] == 'root':
            hit = n

    # sanity check.
    if hit is None:
        raise StandardError

    # return the basal ancestor.
    return hit


def lookup_text(g, text):

    # look to see if we have lu.
    if 'text_lu' in g.graph:
        lu = g.graph['text_lu']

    else:
        # build lookup dictionary.
        lu = {'MainType': {}}
        for n in g.nodes():

            # build the lookup.
            lu[g.node[n]['text']] = n
            if n != 'root':
                mainType = g.node[n]['metamaintype']
                if mainType not in lu['MainType']:
                    lu['MainType'][mainType] = list()
                lu['MainType'][mainType].append(g.node[n]['text'])

        # save it for re-use.
        g.graph['text_lu'] = lu

    # return result. query Main Type first.
    if text in lu['MainType']:
        return lu['MainType'][text]
    elif text in lu:
        return lu[text]
    else:
        return None
