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


		self.node_attributes  = self._get_node_attributes()
		self.group_attributes = self._get_group_attributes()

		self.node_attributes_df = pd.DataFrame(self.node_attributes).set_index('id')

		self.edge_attributes  = self._get_edge_attributes_from_reactions_and_relations()


	######################################
	### Parse tree entries (aka nodes) ###
	######################################
	def _get_node_attributes(self): 

		node_attributes = []

		for entry in self.entries: 

			entry_type = entry.get('type')

			# Parse `gene`-type entries
			if entry_type == 'gene': 
				node_attribute = {
					'type': entry_type, 									# i.e. 'gene'
					'id': int(entry.get('id')), 							# i.e. 84
					'aliases': entry[0].get('name'), 						# i.e. 'ALDOA, ALDA, GSD12, HEL-S-87p...'
					'hsa_tags': entry.get('name'), 							# i.e. 'hsa:226 hsa:229 hsa:230'
					'name': entry[0].get('name').split(', ')[0].rstrip('.') # i.e. 'ALDOA'
				}
				node_attributes.append(node_attribute)
			# Parse `compound`-type entries
			elif entry_type == 'compound': 
				node_attribute = {
					'type': entry_type, 
					'id': int(entry.get('id')),
					'name': entry.get('name')
				}
				node_attributes.append(node_attribute)

			else: pass

		return node_attributes

	def _get_group_attributes(self): 

		group_attributes = []

		for entry in self.entries: 

			if entry.get('type') == 'group':

				group_attribute = {
					'type': entry.get('type'), 
					'id': int(entry.get('id')), 
					'components': [int(x.get('id')) for x in entry.findall('component')]
				}
				group_attributes.append(group_attribute)

		return group_attributes


	######################################
	### Parse tree reactions
	######################################

	def _get_edge_attributes_from_reactions_and_relations(self): 

		# Prioritize reaction edges over relations by populating `edge_attributes` with reaction edges first
		edge_attributes = self._get_edge_attributes_from_reactions()

		existing_edges = [set([edge_attribute['source'], edge_attribute['target']]) for edge_attribute in edge_attributes]
		relation_attributes = self._get_edge_attributes_from_relations()

		for relation_attribute in relation_attributes: 

			if set([relation_attribute['source'], relation_attribute['target']]) not in existing_edges: 

				edge_attributes.append(relation_attribute)

		return edge_attributes


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

		# TODO: Check if edge already exists from reactions function, else return None

		edge_attributes = {
			'source': source,
			'target': target,
			'type': edge_type, 
			'activation': 0,  		# 1 if node1 activates node2, -1 if node2 activates node1, 0 if unknown (ie. oriented==0)
			'oriented': 0,  		# 1 if orientation of arrow is known, 0 otherwise. This may be extraneuous
			'phosphorylation': 0,  
			'glycosylation': 0, 
			'ubiquitination': 0, 
			'methylation': 0
		}

		for interaction in interactions: 
			if interaction in ['activation', 'expression', 'indirect effect']:  edge_attributes.update({ 'activation': 1, 'oriented': 1 })
			elif interaction in ['inhibition', 'repression']:  					edge_attributes.update({ 'activation': -1, 'oriented': 1 })
			elif interaction in ['binding/association']: 						edge_attributes.update({ 'activation': 0, 'oriented': 0 })

			elif interaction == 'phosphorylation': 								edge_attributes.update({ 'phosphorylation': 1, 'oriented': 1 })
			elif interaction == 'dephosphorylation': 							edge_attributes.update({ 'phosphorylation': -1, 'oriented': 1 })				
			elif interaction == 'glycosylation': 								edge_attributes.update({ 'glycosylation': 1, 'oriented': 1 })
			elif interaction == 'ubiquitination': 								edge_attributes.update({ 'ubiquitination': 1, 'oriented': 1 })
			elif interaction == 'methylation': 									edge_attributes.update({ 'methylation': 1, 'oriented': 1 })

			elif interaction == 'compound':										pass
			else: 																pass

		return edge_attributes


	# def _get_edge_attributes_


	# def _get_edge_attributes_from_relations(self): 

	# 	relation_attributes = []

	# 	for relation in self.relations: 

	# 		if relation.get('type') in [ 'PPrel', 'GErel' ]: # TODO: include `PCrel`

	# 			edge_attributes = { 
	# 				'node1': int(relation.get('entry1')),
	# 				'node2': int(relation.get('entry2')),
	# 				'type': relation.get('type'),  # ECrel, PPrel, GErel, PCrel, maplink
	# 				'activation': 0,  # 1 if node1 activates node2, -1 if node2 activates node1, 0 if unknown (ie. oriented==0)
	# 				'oriented': 1,  # 1 if orientation of arrow is known, 0 otherwise. This may be extraneuous
	# 				'phosphorylation': 0,  
	# 				'glycosylation': 0, 
	# 				'ubiquitination': 0, 
	# 				'methylation': 0
	# 			}

	# 			for relation_tag in relation: 

	# 				interaction = relation_tag.get('name')

	# 				if interaction in ['activation', 'expression', 'indirect effect']: 
	# 					edge_attributes['activation'] = 1
	# 				elif interaction in ['inhibition', 'repression']: 
	# 					edge_attributes['activation'] = -1
	# 				elif interaction in ['binding/association']: 
	# 					edge_attributes['activation'] = 0
	# 					edge_attributes['oriented'] = 0

	# 				elif interaction == 'phosphorylation': 
	# 					edge_attributes['phosphorylation'] = 1
	# 				elif interaction == 'dephosphorylation': 
	# 					edge_attributes['phosphorylation'] = -1
						
	# 				elif interaction == 'glycosylation':
	# 					edge_attributes['glycosylation'] = 1
						
	# 				elif interaction == 'ubiquitination': 
	# 					edge_attributes['ubiquitination'] = 1
						
	# 				elif interaction == 'methylation': 
	# 					edge_attributes['methylation'] = 1

	# 			relation_attributes.append(edge_attributes)


	# 				# Force edge through compound node. Unfortunately, these interactions are poorly annotated in KGML.
	# 				# These should be of type 'ECrel', but may be mistakenly annotated as 'PPrel' in KGMLs.
	# 				elif interaction == 'compound': 
	# 					edge1 = {'node1': edge_attributes['node1'], 'node2': int(element.get('value')), 'type': 'PPrel'}
	# 					edge2 = {'node2': edge_attributes['node2'], 'node2': int(element.get('value')), 'type': 'PPrel'}

	# 				else: 
	# 					pass

	# 		# TODO: include ECrel




	# def function to add edges or nodes with optional inputs

	def __str__(self): 

		print_statement = 'Number of entries: {}\nNumber of relations: {}\nNumber of reactions: {}'.format(len(self.entries), len(self.relations), len(self.reactions))
		
		return print_statement