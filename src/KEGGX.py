#!/usr/bin/env python3

import pandas as pd
import networkx as nx

import xml.etree.ElementTree as ET

from itertools import combinations


class KEGGX:

	def __init__(self, KGML_file):

		# Set pathway metadata attributes
		self.root   = ET.parse(KGML_file).getroot()
		self.name   = self.root.get('name') 
		self.org    = self.root.get('org') 
		self.number = self.root.get('number')
		self.title  = self.root.get('title') 
		self.link   = self.root.get('link')

		# Each of the 3 types of elements allowed in KGML files
		self.entries   = self.root.findall('entry')
		self.reactions = self.root.findall('reaction')
		self.relations = self.root.findall('relation')


		self.group_attributes   = self._get_group_attributes()
		self.node_attributes_df = self._get_node_attributes_as_dataframe()
		self.edge_attributes_df = self._get_edge_attributes_as_dataframe()


	######################################
	### Parse tree entries (aka nodes) ###
	######################################
	def _get_node_attributes_as_dataframe(self): 

		node_attributes = []

		for entry in self.entries: 

			entry_type = entry.get('type')
			graphics_attribute = entry.find('graphics').attrib

			# Parse `gene`-type entries
			if entry_type == 'gene': 
				node_attribute = { 'node_type': entry_type, 'id': int(entry.get('id')), 'name': entry[0].get('name').split(', ')[0].rstrip('.') }
				node_attributes.append(node_attribute)
			# Parse `compound`-type entries
			elif entry_type == 'compound': 
				node_attribute = { 'node_type': entry_type, 'id': int(entry.get('id')), 'name': entry.get('name') }
				node_attributes.append(node_attribute)

			else: pass

		node_attributes_df = pd.DataFrame(node_attributes)

		return node_attributes_df


	def _get_group_attributes(self): 

		group_attributes = []

		for entry in self.entries: 

			if entry.get('type') == 'group':

				group_attribute = {
					'node_type': entry.get('type'), 
					'id': int(entry.get('id')), 
					'components': [int(x.get('id')) for x in entry.findall('component')]
				}
				group_attributes.append(group_attribute)

		return group_attributes


	######################################
	### Parse edges
	######################################

	def _get_edge_attributes_as_dataframe(self): 

		# Prioritize reaction edges over relations by populating `edge_attributes` with reaction edges first
		edge_attributes     = self._get_edge_attributes_from_reactions()
		relation_attributes = self._get_edge_attributes_from_relations()

		for relation_attribute in relation_attributes: 

			# Get edges in `edge_attributes` as a list of sets
			existing_edges = [set([edge_attribute['source'], edge_attribute['target']]) for edge_attribute in edge_attributes]

			# Add edge attribute if the edge has not been seen before
			if set([relation_attribute['source'], relation_attribute['target']]) not in existing_edges: 

				edge_attributes.append(relation_attribute)

		# Convert to DataFrame and replace group edges 
		edge_attributes_df = self._replace_group_edges(pd.DataFrame(edge_attributes))

		return edge_attributes_df


	def _get_edge_attributes_from_reactions(self): 

		reaction_attributes = []

		for reaction in self.reactions: 

			compound_id, reaction_name, reaction_type = int(reaction.get('id')), reaction.get('name'), reaction.get('type')
			substrate_ids = [int(substrate.get('id')) for substrate in reaction.findall('substrate')]
			products_ids = [int(product.get('id')) for product in reaction.findall('product')]

			# Add substrate-compound interactions first 
			for substrate_id in substrate_ids: 

				if reaction_type == 'irreversible': 
					reaction_attributes.append(self._populate_edge_attribute(substrate_id, compound_id, reaction_name, ['activation']))
				else: # if `reversible`
					reaction_attributes.append(self._populate_edge_attribute(compound_id, substrate_id, reaction_name, ['activation']))

			# Add compound-product interactions next
			for products_id in products_ids: 
				
				reaction_attributes.append(self._populate_edge_attribute(compound_id, products_id, reaction_name, ['activation']))

		return reaction_attributes


	def _get_edge_attributes_from_relations(self): 

		relation_attributes = []

		for relation in self.relations: 

			source, target = int(relation.get('entry1')), int(relation.get('entry2'))
			edge_type = relation.get('type')
			edge_descriptors = [subtype.get('name') for subtype in relation.findall('subtype')]

			if edge_type in [ 'ECrel', 'PPrel', 'GErel', 'PCrel' ]: 

				# If not compound
				if 'compound' not in edge_descriptors: 

					relation_attributes.append(self._populate_edge_attribute(source, target, edge_type, edge_descriptors))

				else: 
					# Loop through all subtype values, find element that matches the subtype `compound`
					compound_id = [int(subtype.get('value')) for subtype in relation.findall('subtype') if subtype.get('name') == 'compound'][0]

					relation_attributes.append(self._populate_edge_attribute(source, compound_id, edge_type, edge_descriptors))
					relation_attributes.append(self._populate_edge_attribute(compound_id, target, edge_type, edge_descriptors))

		return relation_attributes


	def _populate_edge_attribute(self, source, target, edge_type, interactions):

		# 1 if source activates target, -1 if target activates source, 0 if unknown (ie. oriented==0)
		# 1 if orientation of arrow is known, 0 otherwise. This may be extraneuous
		edge_attributes = { 'source': source, 'target': target, 'type': edge_type, 'activation': 0, 'oriented': 0, 'indirect': 0, 'modification': "" }

		for interaction in interactions: 
			if   interaction in ['activation', 'expression']:	edge_attributes.update({ 'activation':  1, 'oriented': 1 })
			elif interaction in ['inhibition', 'repression']:	edge_attributes.update({ 'activation': -1, 'oriented': 1 })
			elif interaction in ['binding/association']: 		edge_attributes.update({ 'activation':  0, 'oriented': 0 })

			elif interaction == 'indirect effect': 				edge_attributes.update({ 'indirect': 1 })

			elif interaction == 'phosphorylation': 				edge_attributes.update({ 'modification': "+p", 'oriented': 1 })
			elif interaction == 'dephosphorylation': 			edge_attributes.update({ 'modification': "-p", 'oriented': 1 })				
			elif interaction == 'glycosylation': 				edge_attributes.update({ 'modification': "+g", 'oriented': 1 })
			elif interaction == 'ubiquitination': 				edge_attributes.update({ 'modification': "+u", 'oriented': 1 })
			elif interaction == 'methylation': 					edge_attributes.update({ 'modification': "+m", 'oriented': 1 })

			elif interaction == 'compound':						pass
			elif interaction == 'complex': 						edge_attributes.update({ 'activation': 1, 'oriented': 0 }) # complex proteins activate each other bidirectionally
			else: 												pass

		# Add rules here

		return edge_attributes


	def _replace_group_edges(self, edge_attributes_df):

		for group_obj in self.group_attributes: 

			group_id, group_members = group_obj['id'], group_obj['components']
			n_members, n_edges = len(group_members), len(edge_attributes_df)

			# Add edges where `node1` or `node2` is a group member
			for node_type in ['source', 'target']: 
				edges_with_df    = edge_attributes_df[edge_attributes_df[node_type] == group_id]
				edges_without_df = edge_attributes_df[edge_attributes_df[node_type] != group_id]
				# Duplicate rows where `node1` contains the `group_id`
				expanded_edges_df = pd.concat([edges_with_df]*n_members).sort_index() 
				# Replace `node` column with repeating list of `group_members`
				expanded_edges_df[node_type] = group_members*len(edges_with_df)
				# Concatenate the new dataframes and reset index
				edge_attributes_df = pd.concat([expanded_edges_df, edges_without_df]).reset_index(drop=True)
		
			# Add edges *between* group members. 
			group_rows = [ self._populate_edge_attribute(a, b, 'PComplex', ['complex']) for a,b in combinations(group_members, 2) ]
			edge_attributes_df = edge_attributes_df.append(pd.DataFrame(group_rows), ignore_index=True).fillna(0)
		
		return edge_attributes_df
			# print('{} members in group {}. {} edges added.'.format(n_members, group_id, len(self.edge_attributes_df)-n_edges))


	#####
	## Outputs
	#####

	def output_KGML_as_graphml(self, path): 

		# Initialize graph from `edge_attributes_df`
		graph = nx.from_pandas_edgelist(self.edge_attributes_df, 'source', 'target', edge_attr=True, create_using=nx.DiGraph())

		# Retrieve graphics attributes (ie. position, shape, color, etc.)
		node_graphics_attributes = {}

		for entry in self.entries: 
			# Fetch and rename graphics attributes
			node_graphics_attribute = entry.find('graphics').attrib
			node_graphics_attribute['aliases'] = node_graphics_attribute.pop('name', None)
			node_graphics_attribute['shape']   = node_graphics_attribute.pop('type', None)

			node_graphics_attributes[int(entry.get('id'))] = node_graphics_attribute

		# Set node graphics attributes
		nx.set_node_attributes(graph, node_graphics_attributes)
		# Set rest of the node attributes
		nx.set_node_attributes(graph, self.node_attributes_df.set_index('id').to_dict('index'))
		# Relabel so nodes are keyed by gene and compound symbols
		nx.relabel_nodes(graph, {node_id: graph.node[node_id]['name'] for node_id in graph.nodes()}, copy=False)

		nx.write_graphml(graph, path)

		return graph


	def output_KGML_as_networkx(self, directed=False): 

		# Edges with unknown orientation must be reversed and appended to the edge list
		if directed: 
			reverse_edges_df = self.edge_attributes_df[self.edge_attributes_df['oriented'] == 0].rename(columns={'source': 'target', 'target': 'source'})
			bidirected_edge_attributes_df = pd.concat([self.edge_attributes_df, reverse_edges_df])

			graph = nx.from_pandas_edgelist(bidirected_edge_attributes_df, 'source', 'target', edge_attr=True, create_using=nx.DiGraph())

		else: 
			graph = nx.from_pandas_edgelist(self.edge_attributes_df, 'source', 'target', edge_attr=True)

		nx.set_node_attributes(graph, self.node_attributes_df.set_index('id').to_dict('index'))

		nx.relabel_nodes(graph, {node_id: graph.node[node_id]['name'] for node_id in graph.nodes()}, copy=False)

		return graph


	def __str__(self): 

		print_statement = 'Number of entries: {}\nNumber of relations: {}\nNumber of reactions: {}'.format(len(self.entries), len(self.relations), len(self.reactions))
		
		return print_statement