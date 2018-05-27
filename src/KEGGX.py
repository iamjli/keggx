#!/usr/bin/env python3

import pandas as pd
import networkx as nx

import xml.etree.ElementTree as ET

from itertools import combinations, product


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
		self._entries   = self.root.findall('entry')
		self._reactions = self.root.findall('reaction')
		self._relations = self.root.findall('relation')
		self._groups    = self.root.findall('.//entry[@type="group"]')


		self.entry_attributes_df = self._get_entry_attributes_as_dataframe()
		self.node_attributes_df  = self._get_node_attributes_as_dataframe()
		self.edge_attributes_df  = self._get_edge_attributes_as_dataframe()

		self.gene_ids = self._get_node_attributes_as_dataframe(['gene']).index


	#### NODES ####

	def _get_entry_attributes_as_dataframe(self):
		"""
		Parse entry elements and store attributes in a dataframe. Entries include all KEGG pathway 
		elements, incuding maplinks, orthologs, etc.

		Returns: 
			pandas.DataFrame: entry attributes
		"""

		entry_type_df     = pd.DataFrame([entry.attrib for entry in self._entries]).drop(columns=['name', 'link'])
		entry_graphics_df = pd.DataFrame([entry.find('graphics').attrib for entry in self._entries]).rename(columns={'name': 'aliases', 'type': 'shape'})

		entry_attributes_df = pd.concat([entry_type_df, entry_graphics_df], axis=1).fillna("")
		entry_attributes_df['name'] = entry_attributes_df['aliases'].apply(lambda x: x.split(', ')[0].rstrip('.'))
		entry_attributes_df = entry_attributes_df[['id', 'name', 'aliases', 'type', 'x', 'y', 'height', 'width', 'shape', 'bgcolor', 'fgcolor']].set_index('id')

		return entry_attributes_df


	def _get_node_attributes_as_dataframe(self, types=['gene', 'compound']): 
		"""
		Gets entry elements of a particular type.

		Arguments:
			types (list): element types to keep

		Returns: 
			pandas.DataFrame
		"""

		return self.entry_attributes_df[self.entry_attributes_df['type'].isin(types)]


	#### EDGES ####

	def _get_edge_attributes_as_dataframe(self): 
		"""
		Logic for converting reaction and relation elements into a dataframe of edge attributes

		Returns: 
			pandas.DataFrame: edge attributes 
		"""

		edge_attributes_list     = self._get_edge_attributes_from_reactions()
		relation_attributes_list = self._get_edge_attributes_from_relations()

		# Prioritize reaction edges over relations by populating `edge_attributes` with reaction edges first
		for relation_attributes in relation_attributes_list: 

			# Get edges in `edge_attributes` as a list of sets
			existing_edges = [set([edge_attributes['source'], edge_attributes['target']]) for edge_attributes in edge_attributes_list]

			# Add edge attribute if the edge has not been seen before
			if set([relation_attributes['source'], relation_attributes['target']]) not in existing_edges: 

				edge_attributes_list.append(relation_attributes)

		# Convert to DataFrame and replace group edges 
		edge_attributes_df = pd.DataFrame(edge_attributes_list)
		edge_attributes_df = self._replace_group_edges(edge_attributes_df)
		if len(edge_attributes_df) > 0: edge_attributes_df = edge_attributes_df[['source', 'target', 'effect', 'indirect', 'modification', 'type']]

		return edge_attributes_df


	def _get_directed_edge_attributes_as_dataframe(self, edge_attributes_df): 
		# If A<-->B, this function splits into two relations: A-->B and B-->A
		reverse_edges_df = edge_attributes_df[edge_attributes_df['effect'].isin([-2,0,2])].rename(columns={ 'source': 'target', 'target': 'source'})
		directed_edge_attributes_df = pd.concat([edge_attributes_df, reverse_edges_df])

		return directed_edge_attributes_df


	def _get_gene_edge_attributes_as_dataframe(self): 

		compound_ids = self.node_attributes_df[self.node_attributes_df['type'] == 'compound'].index

		directed_edge_attributes_df = self._get_directed_edge_attributes_as_dataframe(self.edge_attributes_df)

		inferred_edges = []
		oriented_edge_attributes_df = directed_edge_attributes_df[directed_edge_attributes_df['effect'] != 0]

		for compound_id in compound_ids: 
			sourced_compounds_df  = oriented_edge_attributes_df[oriented_edge_attributes_df['source'] == compound_id]
			targeted_compounds_df = oriented_edge_attributes_df[oriented_edge_attributes_df['target'] == compound_id]

			source_nodes = targeted_compounds_df['source'].tolist()
			target_nodes = sourced_compounds_df['target'].tolist()

			for source, target in product(source_nodes, target_nodes): 

				if (source in self.gene_ids) & (target in self.gene_ids):

					inferred_edges.append(self._populate_edge_attributes(source, target, "inferred_rxn", ['activation']))

		inferred_edges_df = pd.DataFrame(inferred_edges).drop_duplicates()

		# Remove duplicated edges, consolidate bidirectional edges
		edgelist_as_sets = [set(pair) for pair in inferred_edges_df[['source', 'target']].values]
		inferred_edges_df['effect'] = [edgelist_as_sets.count(pair) for pair in edgelist_as_sets]
		inferred_edges_df = inferred_edges_df[[False if pair in edgelist_as_sets[:i] else True for i,pair in enumerate(edgelist_as_sets)]]

		return inferred_edges_df


	def _populate_edge_attributes(self, source, target, edge_type, interactions): 

		# Attribute `effect` takes values 0 (---), 1 (-->), 2 (<->), or -1 (--|) to indicate cases where
		# orientation is unknown, the edge is activating, the edge is bidirectional (protein complex), or the edge is inhibitory.
		# Perhaps add `binding` as an attribute? Interactions? 
		edge_attributes = { 'source': source, 'target': target, 'type': edge_type, 
							'effect': 0, 'indirect': 0, 'modification': "" }

		# Attributes must be updated in two steps, since descriptors examined in the second loop
		# are more specific than those in the first, and should be used to overwrite them. 
		for interaction in interactions: 

			if   interaction == 'binding/association': edge_attributes.update({ 'effect': 2 })
			elif interaction == 'protein complex':	   edge_attributes.update({ 'effect': 2 }) # not standard type, but including for clarity
			elif interaction == 'bidirected':	   	   edge_attributes.update({ 'effect': 2 }) # not standard type, but including for clarity
			elif interaction == 'dissociation': 	   edge_attributes.update({ 'effect': 1 })
			elif interaction == 'missing interaction': edge_attributes.update({ 'effect': 0 })
			elif interaction == 'indirect effect':     edge_attributes.update({ 'effect': 1, 'indirect': 1 })
			else: pass

		for interaction in interactions: 

			if   interaction == 'phosphorylation':	 edge_attributes.update({ 'effect': 1, 'modification': "+p" })
			elif interaction == 'dephosphorylation': edge_attributes.update({ 'effect': 1, 'modification': "-p" })
			elif interaction == 'glycosylation': 	 edge_attributes.update({ 'effect': 1, 'modification': "+g" })
			elif interaction == 'ubiquitination': 	 edge_attributes.update({ 'effect': 1, 'modification': "+u" })
			elif interaction == 'methylation': 		 edge_attributes.update({ 'effect': 1, 'modification': "+m" })

		for interaction in interactions: 

			if   interaction == 'activation': 	   edge_attributes.update({ 'effect':  1 })
			elif interaction == 'inhibition': 	   edge_attributes.update({ 'effect': -1 })
			elif interaction == 'expression': 	   edge_attributes.update({ 'effect':  1, 'modification': 'e'})
			elif interaction == 'repression': 	   edge_attributes.update({ 'effect': -1, 'modification': 'e'})
			else: pass

		return edge_attributes


	def _get_edge_attributes_from_reactions(self): 

		reaction_attributes_list = []

		for reaction in self._reactions: 

			compound_id, reaction_name, reaction_type = reaction.get('id'), reaction.get('name'), reaction.get('type')
			substrate_ids = [substrate.get('id') for substrate in reaction.findall('substrate')]
			product_ids   = [product.get('id')   for product   in reaction.findall('product')]

			# Add substrate-compound interactions first 
			for substrate_id in substrate_ids: 

				if reaction_type == 'irreversible': 
					reaction_attributes_list.append(self._populate_edge_attributes(substrate_id, compound_id, reaction_name, ['activation']))
				else: 
					reaction_attributes_list.append(self._populate_edge_attributes(compound_id, substrate_id, reaction_name, ['activation']))

			# Add compound-product interactions next
			for product_id in product_ids: 
				
				reaction_attributes_list.append(self._populate_edge_attributes(compound_id, product_id, reaction_name, ['activation']))

		return reaction_attributes_list


	def _get_edge_attributes_from_relations(self): 

		relation_attributes_list = []

		for relation in self._relations: 

			source, target = relation.get('entry1'), relation.get('entry2')
			edge_type = relation.get('type')
			edge_descriptors = [subtype.get('name') for subtype in relation.findall('subtype')]

			# TODO: add support for maplinks?
			if edge_type in [ 'ECrel', 'PPrel', 'GErel', 'PCrel' ]: 

				if 'compound' not in edge_descriptors: 

					relation_attributes_list.append(self._populate_edge_attributes(source, target, edge_type, edge_descriptors))

				else: 
					# Get compound id by first searching through subtypes with `name` attribute equal to 'compound', then retrieving 'value'
					compound_id = relation.find('.//subtype[@name="compound"]').get('value')

					relation_attributes_list.append(self._populate_edge_attributes(source, compound_id, edge_type, edge_descriptors))
					relation_attributes_list.append(self._populate_edge_attributes(compound_id, target, edge_type, edge_descriptors))

		return relation_attributes_list


	def _replace_group_edges(self, edge_attributes_df):

		for group_element in self._groups: 

			group_id = group_element.get('id')
			group_members = [component.get('id') for component in group_element.findall('component')]

			# TODO: Make sure these edges haven't been added yet.
			# Add edges where `node1` or `node2` is a group member
			for node_type in ['source', 'target']: 
				edges_with_df    = edge_attributes_df[edge_attributes_df[node_type] == group_id]
				edges_without_df = edge_attributes_df[edge_attributes_df[node_type] != group_id]
				# Duplicate rows where `node1` contains the `group_id`
				expanded_edges_df = pd.concat([edges_with_df]*len(group_members)).sort_index() 
				# Replace `node` column with repeating list of `group_members`
				expanded_edges_df[node_type] = group_members*len(edges_with_df)
				# Concatenate the new dataframes and reset index
				edge_attributes_df = pd.concat([expanded_edges_df, edges_without_df]).reset_index(drop=True)
		
			# Add edges *between* group members. Complexed proteins are essentially `binding/association`
			group_rows = [ self._populate_edge_attributes(a, b, 'PComplex', ['protein complex']) for a,b in combinations(group_members, 2) ]
			edge_attributes_df = edge_attributes_df.append(pd.DataFrame(group_rows), ignore_index=True).fillna(0)
		
		return edge_attributes_df


	#### OUTPUTS ####

	def output_KGML_as_graphml(self, path, display='full'): 

		# Initialize graph from `edge_attributes_df`
		if len(self.edge_attributes_df) > 0:
			graph = nx.from_pandas_edgelist(self.edge_attributes_df, 'source', 'target', edge_attr=True, create_using=nx.DiGraph())
		else: 
			graph = nx.DiGraph()

		# Detailed visualization includes singletons as non-gene or compound nodes, such as orthology, titles, etc.
		if display == 'full': 
			graph.add_nodes_from(self.entry_attributes_df.index)
		elif display == 'genes_and_compounds': 
			pass
		elif display == 'genes': 
			graph = nx.DiGraph(graph.subgraph(self.node_attributes_df.index[self.node_attributes_df['type'] == 'gene']))

			inferred_edges_df = self._get_gene_edge_attributes_as_dataframe()
			inferred_edges_graph = nx.from_pandas_edgelist(inferred_edges_df, 'source', 'target', edge_attr=True, create_using=nx.DiGraph())

			graph = nx.compose(graph, inferred_edges_graph)

		nx.set_node_attributes(graph, self.entry_attributes_df.to_dict('index'))

		nx.write_graphml(graph, path)

		return graph
		

	def output_KGML_as_directed_networkx(self, genes_only=True): 

		directed_edge_attributes_df = self._get_directed_edge_attributes_as_dataframe(self.edge_attributes_df)

		if genes_only: 

			inferred_edge_attributes_df = self._get_directed_edge_attributes_as_dataframe(self._get_gene_edge_attributes_as_dataframe())
			directed_edge_attributes_df = pd.concat([directed_edge_attributes_df, inferred_edge_attributes_df])

			graph = nx.from_pandas_edgelist(directed_edge_attributes_df, 'source', 'target', edge_attr=True, create_using=nx.DiGraph())
			graph = nx.DiGraph(graph.subgraph(self.node_attributes_df.index[self.node_attributes_df['type'] == 'gene']))

		else: 
			
			graph = nx.from_pandas_edgelist(directed_edge_attributes_df, 'source', 'target', edge_attr=True, create_using=nx.DiGraph())

		nx.set_node_attributes(graph, self.entry_attributes_df.to_dict('index'))
		nx.relabel_nodes(graph, { node_id: graph.node[node_id]['name'] for node_id in graph.nodes() }, copy=False)

		return graph
