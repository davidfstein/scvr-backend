import scanpy as sc
from scvr import converters
from scipy.sparse import isspmatrix
from pandas.api.types import is_string_dtype,is_numeric_dtype
from glob import glob
import pathlib
import os
from flask import Flask, jsonify, request
from flask_cors import CORS, cross_origin
import matplotlib as mpl
import networkx as nx
import numpy as np
#import stream as st
import pandas as pd
import fetch_data
import resource

# Initialize app
APP_PATH = str(pathlib.Path(__file__).parent.resolve())
DATASET_DIRECTORY = os.path.join(APP_PATH, "app_datasets")
UPLOAD_DIRECTORY = os.path.join(APP_PATH, "app_uploaded_files")
QR_DIRECTORY = os.path.join(APP_PATH, "assets")
server = Flask(__name__)
cors = CORS(server)
server.config['CORS_HEADERS'] = 'Content-Type'

if not os.path.exists(DATASET_DIRECTORY):
    os.makedirs(DATASET_DIRECTORY)

# Quick and dirty
def get_tool_type(file):
    for tool in ['stream', 'paga', 'scanpy', 'seurat']:
        if tool in file.lower():
            return tool

@server.route("/databases", methods=["GET", "POST"])
def get_databases():
    adata_files = glob(os.path.join(DATASET_DIRECTORY, "*"))
    adata_list = list()
    for file in adata_files:
        tool = get_tool_type(file)
        label = os.path.basename(file).split(".")[0]
        label_list = label.split("_")
        adata_list.append(
            {
                "value": label_list[0],
                "label": label,
                "path": file,
                "type": tool,
            }
        )
    return jsonify(adata_list)


@server.route("/data_type", methods=["GET", "POST"])
def get_dataset_type():
    """
    http://127.0.0.1:8000/data_type?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    return jsonify({"type": db_name.split("_")[1]})


def get_dataset_type_adata(db_name):
    return db_name.split("_")[1]


@server.route("/coordinates", methods=["GET", "POST"])
def get_coordinates():
    """
    http://127.0.0.1:8000/coordinates?db_name=1_scanpy_10xpbmc&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=3_velocity_pancrease&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=4_seurat_10xpbmc&embed=umap
    http://127.0.0.1:8000/coordinates?db_name=5_stream_nestorowa16&embed=umap
    """
    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    if get_dataset_type_adata(db_name).lower() in ["scanpy", "velocity", "seurat", "paga"]:
        adata = sc.read(filename)
        embed = request.args.get("embed")
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")

    list_cells = []
    for i in range(adata.shape[0]):
        dict_coord_cells = dict()
        dict_coord_cells["cell_id"] = adata.obs_names[i]
        if get_dataset_type_adata(db_name).lower() in ["scanpy", "paga", "velocity"]:
            dict_coord_cells["x"] = str(adata.obsm[f"X_{embed}"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm[f"X_{embed}"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm[f"X_{embed}"][i, 2])
        elif get_dataset_type_adata(db_name).lower() == "seurat":
            dict_coord_cells["x"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 0])
            dict_coord_cells["y"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 1])
            dict_coord_cells["z"] = str(adata.obsm[f"{embed}_cell_embeddings"][i, 2])
        elif get_dataset_type_adata(db_name).lower() == "stream":
            file_path = os.path.join(adata.uns["workdir"], "test")
            if not os.path.exists(file_path):
                os.makedirs(file_path)
            flat_tree = adata.uns["flat_tree"]
            epg = adata.uns["epg"]
            epg_node_pos = nx.get_node_attributes(epg, "pos")
            ft_node_label = nx.get_node_attributes(flat_tree, "label")
            ft_node_pos = nx.get_node_attributes(flat_tree, "pos")
            list_curves = []
            for edge_i in flat_tree.edges():
                branch_i_pos = np.array(
                    [epg_node_pos[i] for i in flat_tree.edges[edge_i]["nodes"]]
                )
                df_coord_curve_i = pd.DataFrame(branch_i_pos)
                dict_coord_curves = dict()
                dict_coord_curves["branch_id"] = (
                    ft_node_label[edge_i[0]] + "_" + ft_node_label[edge_i[1]]
                )
                dict_coord_curves["xyz"] = [
                    {
                        "x": df_coord_curve_i.iloc[j, 0],
                        "y": df_coord_curve_i.iloc[j, 1],
                        "z": df_coord_curve_i.iloc[j, 2],
                    }
                    for j in range(df_coord_curve_i.shape[0])
                ]
                list_curves.append(dict_coord_curves)

            ## output topology of stream graph
            dict_nodes = dict()
            list_edges = []
            for node_i in flat_tree.nodes():
                dict_nodes_i = dict()
                dict_nodes_i["node_name"] = ft_node_label[node_i]
                dict_nodes_i["xyz"] = {
                    "x": ft_node_pos[node_i][0],
                    "y": ft_node_pos[node_i][1],
                    "z": ft_node_pos[node_i][2],
                }
                dict_nodes[ft_node_label[node_i]] = dict_nodes_i
            for edge_i in flat_tree.edges():
                dict_edges = dict()
                dict_edges["nodes"] = [
                    ft_node_label[edge_i[0]],
                    ft_node_label[edge_i[1]],
                ]
                dict_edges["weight"] = 1
                list_edges.append(dict_edges)

            list_cells = []
            for i in range(adata.shape[0]):
                dict_coord_cells = dict()
                dict_coord_cells['cell_id'] = adata.obs_names[i]
                dict_coord_cells['x'] = adata.obsm['X_dr'][i,0]
                dict_coord_cells['y'] = adata.obsm['X_dr'][i,1]
                dict_coord_cells['z'] = adata.obsm['X_dr'][i,2]
                list_cells.append(dict_coord_cells)
            return jsonify(
                {"nodes": dict_nodes, "edges": list_edges, "graph": list_curves, "cells": list_cells}
            )
        else:
            raise TypeError("not supported format")
        list_cells.append(dict_coord_cells)
    return jsonify(list_cells)


@server.route("/features", methods=["GET", "POST"])
def get_features():
    """
    scanpy examples:
      http://127.0.0.1:8000/features?db_name=1_scanpy_10xpbmc&feature=louvain
      http://127.0.0.1:8000/features?db_name=1_scanpy_10xpbmc&feature=expression&gene=SUMO3

    seurat examples:
      http://127.0.0.1:8000/features?db_name=4_seurat_10xpbmc&feature=expression&gene=SUMO3
      http://127.0.0.1:8000/features?db_name=4_seurat_10xpbmc&feature=expression&gene=SUMO3

    velocity examples:
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=clusters
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=expression&gene=Rbbp7
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity&embed=umap&time=None
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity&embed=umap&time=1
      http://127.0.0.1:8000/features?db_name=3_velocity_pancrease&feature=velocity&embed=umap&time=10
    """
    database = request.args.get("db_name")
    feature = request.args.get("feature")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{database}.*"))[0]

    db_type = get_dataset_type_adata(filename)
    if feature.lower() == "velocity":
        embed = request.args.get("embed")

    if get_dataset_type_adata(database).lower() in ["scanpy", "velocity", "seurat", "paga"]:
        adata = sc.read(filename)
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")

    list_metadata = []
    if feature in get_available_annotations_adata(adata):  # cluster columns
        adata.obs[feature] = adata.obs[feature].astype('category')
        if f"{feature}_colors" in adata.uns.keys():
            dict_colors = {
                feature: dict(
                    zip(adata.obs[feature].cat.categories, adata.uns[f"{feature}_colors"])
                )
            }
        else:
            dict_colors = {
                feature: dict(
                    zip(adata.obs[feature].cat.categories, get_colors(adata, feature))
                )
            }
        for i in range(adata.shape[0]):
            dict_metadata = dict()
            dict_metadata["cell_id"] = adata.obs_names[i]
            dict_metadata["label"] = adata.obs[feature].tolist()[i]
            dict_metadata["clusters"] = adata.obs[feature].tolist()[i]
            dict_metadata["clusters_color"] = dict_colors[feature][
                dict_metadata["clusters"]
            ]
            list_metadata.append(dict_metadata)
    elif feature in ["expression", "rna"]:  # pseudotime or latent_time columns
        gene = request.args.get("gene")
        if gene not in adata.var_names:
            return jsonify({})
        else:
            if "time" in feature:
                values = adata.obs[feature]
            else:
                if db_type == "seurat":
                    values = (
                        adata[:, gene].layers["norm_data"].toarray()[:, 0]
                        if isspmatrix(adata.layers["norm_data"])
                        else adata[:, gene].layers["norm_data"][:, 0]
                    )
                else:
                    values = (
                        adata[:, gene].X.toarray()[:, 0]
                        if isspmatrix(adata.X)
                        else adata[:, gene].X[:, 0]
                    )

            cm = mpl.cm.get_cmap("viridis", 512)
            norm = mpl.colors.Normalize(vmin=0, vmax=max(values), clip=True)
            list_metadata = []
            for i, x in enumerate(adata.obs_names):
                dict_genes = dict()
                dict_genes["cell_id"] = x
                dict_genes["color"] = mpl.colors.to_hex(cm(norm(values[i])))
                list_metadata.append(dict_genes)
    elif feature == "velocity":
        list_metadata = []
        time = request.args.get("time")
        for i in range(adata.shape[0]):
            dict_coord_cells = dict()
            if isinstance(adata.obs_names[i], bytes):
                dict_coord_cells["cell_id"] = adata.obs_names[i].decode("utf-8")
            else:
                dict_coord_cells["cell_id"] = adata.obs_names[i]

            dict_coord_cells["x0"] = str(adata.obsm[f"X_{embed}"][i, 0])
            dict_coord_cells["y0"] = str(adata.obsm[f"X_{embed}"][i, 1])
            dict_coord_cells["z0"] = str(adata.obsm[f"X_{embed}"][i, 2])

            if time == "None":
                dict_coord_cells["x1"] = str(adata.obsm[f"velocity_{embed}"][i, 0])
                dict_coord_cells["y1"] = str(adata.obsm[f"velocity_{embed}"][i, 1])
                dict_coord_cells["z1"] = str(adata.obsm[f"velocity_{embed}"][i, 2])
            elif time in list(map(str, [0.01, 0.1, 1, 5, 10, 20, 30, 50, 100])):
                dict_coord_cells["x1"] = str(
                    adata.obsm[f"absolute_velocity_{embed}_{time}s"][i, 0]
                )
                dict_coord_cells["y1"] = str(
                    adata.obsm[f"absolute_velocity_{embed}_{time}s"][i, 1]
                )
                dict_coord_cells["z1"] = str(
                    adata.obsm[f"absolute_velocity_{embed}_{time}s"][i, 2]
                )
            else:
                return jsonify({})
            list_metadata.append(dict_coord_cells)
    elif feature == "paga":
        G = nx.from_numpy_matrix(adata.uns["paga"]["connectivities"].toarray())
        adata.uns["paga"]["pos"] = get_paga3d_pos(adata)
        ## output coordinates of paga graph
        list_lines = []
        for edge_i in G.edges():
            dict_coord_lines = dict()
            dict_coord_lines["branch_id"] = [[str(edge_i[0]), str(edge_i[1])]]
            dict_coord_lines["xyz"] = [
                {"x": pos[0], "y": pos[1], "z": pos[2]}
                for pos in adata.uns["paga"]["pos"][[edge_i[0], edge_i[1]], :]
            ]
            list_lines.append(dict_coord_lines)

        ## output topology of paga graph
        dict_nodes = dict()
        list_edges = []
        dict_nodename = {
            i: adata.obs[adata.uns["paga"]["groups"]].cat.categories[i]
            for i in G.nodes()
        }
        for node_i in G.nodes():
            dict_nodes_i = dict()
            dict_nodes_i["node_name"] = dict_nodename[node_i]
            dict_nodes_i["xyz"] = {
                "x": adata.uns["paga"]["pos"][:, 0][node_i],
                "y": adata.uns["paga"]["pos"][:, 1][node_i],
                "z": adata.uns["paga"]["pos"][:, 2][node_i],
            }
            dict_nodes[node_i] = dict_nodes_i
        for edge_i in G.edges():
            dict_edges = dict()
            dict_edges["nodes"] = [str(edge_i[0]), str(edge_i[1])]
            dict_edges["weight"] = adata.uns["paga"]["connectivities"][
                edge_i[0], edge_i[1]
            ]
            list_edges.append(dict_edges)
        list_metadata = {"nodes": dict_nodes, "edges": list_edges}
    elif feature == "curves":
        flat_tree = adata.uns['flat_tree']
        epg = adata.uns['epg']
        epg_node_pos = nx.get_node_attributes(epg,'pos')
        ft_node_label = nx.get_node_attributes(flat_tree,'label')
        ft_node_pos = nx.get_node_attributes(flat_tree,'pos')
        list_curves = []
        for edge_i in flat_tree.edges():
            branch_i_pos = np.array([epg_node_pos[i] for i in flat_tree.edges[edge_i]['nodes']])
            df_coord_curve_i = pd.DataFrame(branch_i_pos)
            dict_coord_curves = dict()
            dict_coord_curves['branch_id'] = ft_node_label[edge_i[0]] + '_' + ft_node_label[edge_i[1]]
            dict_coord_curves['xyz'] = [{'x':df_coord_curve_i.iloc[j,0],
                                         'y':df_coord_curve_i.iloc[j,1],
                                         'z':df_coord_curve_i.iloc[j,2]} for j in range(df_coord_curve_i.shape[0])]
            list_curves.append(dict_coord_curves)
        list_metadata = list_curves
    
    return jsonify({feature: list_metadata})


def get_paga3d_pos(adata):
    assert (
        adata.obsm["X_umap"].shape[1] >= 3
    ), """The embedding space should have at least three dimensions. 
        please set `n_component = 3` in `sc.tl.umap()`"""
    groups = adata.obs[adata.uns["paga"]["groups"]]
    connectivities_coarse = adata.uns["paga"]["connectivities"]
    paga3d_pos = np.zeros((connectivities_coarse.shape[0], 3))
    for i in range(connectivities_coarse.shape[0]):
        subset = (groups == groups.cat.categories[i]).values
        paga3d_pos[i] = np.median(adata.obsm["X_umap"][subset], axis=0)
    return paga3d_pos


@server.route("/columns", methods=["GET", "POST"])
def get_available_annotations():
    """
    http://127.0.0.1:8000/columns?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    if get_dataset_type_adata(db_name).lower() in ["scanpy", "velocity", "seurat", "paga"]:
        adata = sc.read(filename)
    else:
        adata = st.read(filename, file_format="pkl", workdir="./")
    return jsonify(list(adata.obs.columns))


def get_available_annotations_adata(adata):
    return adata.obs.columns


@server.route("/genes", methods=["GET", "POST"])
def get_genes():
    """
    http://127.0.0.1:8000/genes?db_name=1_scanpy_10xpbmc
    """
    db_name = request.args.get("db_name")
    adata = None
    if get_dataset_type_adata(db_name).lower() == 'stream':
        adata = st.read(glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0], file_format="pkl", workdir="./")
    else:
        adata = sc.read(glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0])
    return jsonify(list(adata.var_names))


def get_genes_adata(adata):
    return adata.var_names


@server.route("/ts", methods=["GET", "POST"])
def get_ts():
    """
    velocity examples:
    http://127.0.0.1:8000/ts?db_name=3_velocity_pancrease&feature=clusters
    """

    db_name = request.args.get("db_name")
    filename = glob(os.path.join(DATASET_DIRECTORY, f"{db_name}.*"))[0]

    if get_dataset_type_adata(db_name) == "velocity":
        adata = sc.read(filename)
        ts = [k.replace('absolute_velocity_umap_', '').replace('s', '')
              for k in adata.obsm.keys() if k.startswith('absolute')]
    return jsonify(list(ts))

def get_colors(adata,ann):
    df_cell_colors = pd.DataFrame(index=adata.obs.index)
    if(is_numeric_dtype(adata.obs[ann])):
        cm = mpl.cm.get_cmap('viridis',512)
        norm = mpl.colors.Normalize(vmin=0, vmax=max(adata.obs[ann]),clip=True)
        df_cell_colors[ann+'_color'] = [mpl.colors.to_hex(cm(norm(x))) for x in adata.obs[ann]]
    else:
        df_cell_colors[ann+'_color'] = ''

        adata.obs[ann] = adata.obs[ann].astype('category')
        categories = adata.obs[ann].cat.categories
        length = len(categories)
        # check if default matplotlib palette has enough colors
        # mpl.style.use('default')
        if len(mpl.rcParams['axes.prop_cycle'].by_key()['color']) >= length:
            cc = mpl.rcParams['axes.prop_cycle']()
            palette = [next(cc)['color'] for _ in range(length)]
        else:
            if length <= 20:
                palette = default_20
            elif length <= 28:
                palette = default_28
            elif length <= len(default_102):  # 103 colors
                palette = default_102
            else:
                rgb_rainbow = mpl.cm.rainbow(np.linspace(0,1,length))
                palette = [mpl.colors.rgb2hex(rgb_rainbow[i,:-1]) for i in range(length)]
        for i,x in enumerate(categories):
            id_cells = np.where(adata.obs[ann]==x)[0]
            df_cell_colors.loc[df_cell_colors.index[id_cells],ann+'_color'] = palette[i]
    return(df_cell_colors[ann+'_color'].tolist())


### modifed from scanpy palettes https://github.com/theislab/scanpy/blob/master/scanpy/plotting/palettes.py

"""Color palettes in addition to matplotlib's palettes."""

from matplotlib import cm, colors

# Colorblindness adjusted vega_10
# See https://github.com/theislab/scanpy/issues/387
vega_10 = list(map(colors.to_hex, cm.tab10.colors))
vega_10_scanpy = vega_10.copy()
vega_10_scanpy[2] = '#279e68'  # green
vega_10_scanpy[4] = '#aa40fc'  # purple
vega_10_scanpy[8] = '#b5bd61'  # kakhi

# default matplotlib 2.0 palette
# see 'category20' on https://github.com/vega/vega/wiki/Scales#scale-range-literals
vega_20 = list(map(colors.to_hex, cm.tab20.colors))

# reorderd, some removed, some added
vega_20_scanpy = [
    *vega_20[0:14:2], *vega_20[16::2],  # dark without grey
    *vega_20[1:15:2], *vega_20[17::2],  # light without grey
    '#ad494a', '#8c6d31',  # manual additions
]
vega_20_scanpy[2] = vega_10_scanpy[2]
vega_20_scanpy[4] = vega_10_scanpy[4]
vega_20_scanpy[7] = vega_10_scanpy[8]  # kakhi shifted by missing grey
# TODO: also replace pale colors if necessary

default_20 = vega_20_scanpy

# https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
# update 1
# orig reference http://epub.wu.ac.at/1692/1/document.pdf
zeileis_28 = [
    "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3",
    "#8595e1", "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593",
    "#c6dec7", "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6", "#d5eae7",
    "#f3e1eb", "#f6c4e1", "#f79cd4",
    '#7f7f7f', "#c7c7c7", "#1CE6FF", "#336600",  # these last ones were added,
]

default_28 = zeileis_28

# from http://godsnotwheregodsnot.blogspot.de/2012/09/color-distribution-methodology.html
godsnot_102 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72",
]

default_102 = godsnot_102
