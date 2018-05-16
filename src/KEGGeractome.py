#!/usr/bin/env python3

import pandas as pd
import networkx as nx

import xml.etree.ElementTree as ET

from itertools import combinations


class KEGGeractome:

	def __init__(self, KGML_file):

		self.root   = ET.parse(KGML_file).getroot()
		self.name   = self.root.get("name") 
		self.org    = self.root.get("org") 
		self.number = self.root.get("number")
		self.title  = self.root.get("title") 
		self.link   = self.root.get("link")

		self.node_attributes    = [self._node_attribute_from_gene_entry(entry) for entry in self.root.findall("entry") if entry.get("type") == "gene"]
		self.node_attributes_df = pd.DataFrame(self.node_attributes).set_index("id")
		self.complex_attributes = [self._complex_attributes_from_group_entry(entry) for entry in self.root.findall("entry") if entry.get("type") == "group"]

		self.edge_attributes    = [self._edge_attributes_from_PP_relation(relation) for relation in self.root.findall("relation") if relation.get("type") in ["PPrel", "PCrel"]]
		self.edge_attributes_df = pd.DataFrame(self.edge_attributes)

		self._replace_group_edges()

		self.graph, self.digraph = self._generate_graph_from_attributes()


	## Parse KGML entries
	def _node_attribute_from_gene_entry(self, entry): 

		node_attribute = {
			"type": entry.get("type"), # i.e. 'gene'
			"id": int(entry.get("id")), # i.e. 84
			"aliases": entry[0].get("name"), # i.e. 'ALDOA, ALDA, GSD12, HEL-S-87p...'
			"hsa_tags": entry.get("name"), # i.e. 'hsa:226 hsa:229 hsa:230'
			"first_name": entry[0].get("name").split(', ')[0].rstrip('.') # i.e. 'ALDOA'
		}

		return node_attribute

	def _complex_attributes_from_group_entry(self, entry): 

		complex_attributes = {
			"type": entry.get("type"), # i.e. 'gene'
			"id": int(entry.get("id")), # i.e. 84
			"components": [int(x.get("id")) for x in entry.findall("component")]
		}

		return complex_attributes

	def _edge_attributes_from_PP_relation(self, relation): 

		# directed, sign, 
		edge_attributes = { 
			"node1": int(relation.get("entry1")),
			"node2": int(relation.get("entry2")),
			"type": relation.get("type"),  # ECrel, PPrel, GErel, PCrel, maplink
			"activation": 0,  # 1 if node1 activates node2, -1 if node2 activates node1, 0 if unknown (ie. oriented==0)
			"oriented": 1,  # 1 if orientation of arrow is known, 0 otherwise. This may be extraneuous
			"phosphorylation": 0,  
			"glycosylation": 0, 
			"ubiquitination": 0, 
			"methylation": 0
		}

		for element in relation: 

			interaction = element.get("name")

			if interaction in ["activation", "expression", "indirect effect"]: 
				edge_attributes["activation"] = 1
			elif interaction in ["inhibition", "repression"]: 
				edge_attributes["activation"] = -1
			elif interaction in ["binding/association"]: 
				edge_attributes["activation"] = 0
				edge_attributes["oriented"] = 0

			elif interaction == "phosphorylation": 
				edge_attributes["phosphorylation"] = 1
			elif interaction == "dephosphorylation": 
				edge_attributes["phosphorylation"] = -1
				
			elif interaction == "glycosylation":
				edge_attributes["glycosylation"] = 1
				
			elif interaction == "ubiquitination": 
				edge_attributes["ubiquitination"] = 1
				
			elif interaction == "methylation": 
				edge_attributes["methylation"] = 1

			# Force edge through compound node. Unfortunately, these interactions are poorly annotated in KGML.
			# These should be of type "ECrel", but may be mistakenly annotated as "PPrel" in KGMLs.
			elif interaction == "compound": 
				edge1 = {"node1": edge_attributes["node1"], "node2": int(element.get("value")), "type": "PPrel"}
				edge2 = {"node2": edge_attributes["node2"], "node2": int(element.get("value")), "type": "PPrel"}

			else: 
				pass

		return edge_attributes


	def _replace_group_edges(self):

		for group_obj in self.complex_attributes: 
			group_id = group_obj["id"] 
			group_members = group_obj["components"]
			n_members, n_edges = len(group_members), len(self.edge_attributes_df)

			# Add edges where `node1` or `node2` is a group member
			for node_type in ["node1", "node2"]: 
				edges_with_df = self.edge_attributes_df[self.edge_attributes_df[node_type] == group_id]
				edges_without_df = self.edge_attributes_df[self.edge_attributes_df[node_type] != group_id]
				# Duplicate rows where `node1` contains the `group_id`
				expanded_edges_df = pd.concat([edges_with_df]*n_members).sort_index() 
				# Replace `node` column with repeating list of `group_members`
				expanded_edges_df[node_type] = group_members*len(edges_with_df)
				# Concatenate the new dataframes and reset index
				self.edge_attributes_df = pd.concat([expanded_edges_df, edges_without_df]).reset_index(drop=True)
		
			# Add edges *between* group members. `sign` will implicity be assigned 0 to indicate no directionality. 
			group_rows = [{"node1": a, "node2": b, "type": "complex", "activation": 0, "oriented": 0} for a,b in combinations(group_members, 2)]
			self.edge_attributes_df = self.edge_attributes_df.append(pd.DataFrame(group_rows), ignore_index=True).fillna(0)
		
			# print("{} members in group {}. {} edges added.".format(n_members, group_id, len(self.edge_attributes_df)-n_edges))


	def _generate_graph_from_attributes(self): 

		# Get node dataframe as a list of dictionaries keyed by node
		nodes_keyed_by_node = [(key, dic) for key, dic in self.node_attributes_df.to_dict("index").items()]
		# Get edge dataframe as a list of dictionaries keyed by nodes. In order to accomodate forward and reverse edges associated 
		# with undirected edges, we filter for edges that are undirected, then append their reverse direction to the list
		edge_attribute_columns = [col for col in self.edge_attributes_df.columns if not col.startswith("node")]
		directed_edges_keyed_by_nodes = [(row["node1"], row["node2"], row[edge_attribute_columns].to_dict()) for _,row in self.edge_attributes_df.iterrows()]
		reversed_edges_keyed_by_nodes = [(row["node2"], row["node1"], row[edge_attribute_columns].to_dict()) for _,row in self.edge_attributes_df.iterrows() if row["oriented"] == 0]
		all_edges_keyed_by_nodes = [*directed_edges_keyed_by_nodes, *reversed_edges_keyed_by_nodes]

		# Create undirected graph. 
		graph = nx.Graph()
		graph.add_nodes_from(nodes_keyed_by_node)
		graph.add_edges_from(directed_edges_keyed_by_nodes)
		nx.relabel_nodes(graph, {node_id: graph.node[node_id]["first_name"] for node_id in graph.nodes()}, copy=False)

		# Create directed graph. 
		digraph = nx.DiGraph()
		digraph.add_nodes_from(nodes_keyed_by_node)
		digraph.add_edges_from(all_edges_keyed_by_nodes)
		nx.relabel_nodes(digraph, {node_id: digraph.node[node_id]["first_name"] for node_id in digraph.nodes()}, copy=False)

		return graph, digraph


	def get_edges_from_protein_list(self, proteins): 

		subgraph = self.graph.subgraph(proteins)
		subgraph_edges = subgraph.edges(data=True)

		# for _,_,attr in subgraph_edges: 
		# 	attr["pathway_id"] = self.number

		return subgraph_edges

	def get_path_between_proteins(self, protein1, protein2, max_len=2): 

		pass





"""
attribute value		explanation
ortholog			the node is a KO (ortholog group)
enzyme				the node is an enzyme
reaction			the node is a reaction
gene				the node is a gene product (mostly a protein)
group				the node is a complex of gene products (mostly a protein complex)
compound			the node is a chemical compound (including a glycan)
map					the node is a linked pathway map
brite				the node is a linked brite hierarchy
other				the node is an unclassified type
"""