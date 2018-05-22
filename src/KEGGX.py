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

		self.edge_attributes  = self._get_edge_attributes_from_reactions()


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
					'type': entry.get('type'), # i.e. 'gene'
					'id': int(entry.get('id')), # i.e. 84
					'components': [int(x.get('id')) for x in entry.findall('component')]
				}
				group_attributes.append(group_attribute)

		return group_attributes


	######################################
	### Parse tree reactions
	######################################
	def _get_edge_attributes_from_reactions(self): 

		reaction_attributes = []

		for reaction in self.reactions: 

			compound_id, reaction_name, reaction_type = int(reaction.get('id')), reaction.get('name'), reaction.get('type')
			substrate_ids = [int(substrate.get('id')) for substrate in reaction.findall('substrate')]
			products_ids = [int(product.get('id')) for product in reaction.findall('product')]

			# Add substrate-compound interactions first 
			for substrate_id in substrate_ids: 

				if reaction_type == 'irreversible': 
					reaction_attribute = {
						'source': substrate_id,
						'target': compound_id,
						'reaction_name': reaction_name, 
						'oriented': 1
					}
				else: # if `reversible`
					reaction_attribute = {
						'source': compound_id,
						'target': substrate_id,
						'reaction_name': reaction_name, 
						'oriented': 1
					}
				reaction_attributes.append(reaction_attribute)

			# Add compound-product interactions next
			for products_id in products_ids: 
				reaction_attribute = {
					'source': compound_id,
					'target': products_id,
					'reaction_name': reaction_name, 
					'oriented': 1
				}
				reaction_attributes.append(reaction_attribute)

		return reaction_attributes

	# def function to add edges or nodes with optional inputs

	def __str__(self): 

		print_statement = "Number of entries: {}\nNumber of relations: {}\nNumber of reactions: {}".format(len(self.entries), len(self.relations), len(self.reactions))
		
		return print_statement