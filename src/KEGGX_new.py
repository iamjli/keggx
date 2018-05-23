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



	#### EDGES ####

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
			elif interaction == 'dissociation': 	   edge_attributes.update({ 'effect': 1 })
			elif interaction == 'missing interaction': edge_attributes.update({ 'effect': 0 })
			else: pass

		for interaction in interactions: 

			if   interaction == 'activation': 	   edge_attributes.update({ 'effect':  1 })
			elif interaction == 'inhibition': 	   edge_attributes.update({ 'effect': -1 })
			elif interaction == 'expression': 	   edge_attributes.update({ 'effect':  1, 'modification': 'e'})
			elif interaction == 'repression': 	   edge_attributes.update({ 'effect': -1, 'modification': 'e'})
			elif interaction == 'indirect effect': edge_attributes.update({ 'indirect': 1 })
			else: pass

		return edge_attributes


	def _get_edge_attributes_from_reaction(self): 

		reaction_attributes_list = []

		for reaction in self.reactions: 

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

		for relation in self.relations: 

			source, target = relation.get('entry1'), relation.get('entry2')
			edge_type = relation.get('type')
			edge_descriptors = [subtype.get('name') for subtype in relation.findall('subtype')]

			# TODO: add support for maplinks?
			if edge_type in [ 'ECrel', 'PPrel', 'GErel', 'PCrel' ]: 

				if 'compound' not in edge_descriptors: 

					relation_attributes_list.append(self._populate_edge_attributes(source, target, edge_type, edge_descriptors))

				else: 

					# Get compound id by first searching through subtypes with `name` attribute equal to 'compound', then retrieving 'value'
					compound_id = relation.find(".//subtype[@name='compound']").get('value')

					relation_attributes_list.append(self._populate_edge_attributes(source, compound_id, edge_type, edge_descriptors))
					relation_attributes_list.append(self._populate_edge_attributes(compound_id, target, edge_type, edge_descriptors))

		return relation_attributes_list


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



