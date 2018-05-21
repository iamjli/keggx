#!/usr/bin/env python3

import pandas as pd
import networkx as nx

import xml.etree.ElementTree as ET

from itertools import combinations


class KEGGX:

	def __init__(self, KGML_file):

		# Set pathway attributes
		self.root   = ET.parse(KGML_file).getroot()
		self.name   = self.root.get("name") 
		self.org    = self.root.get("org") 
		self.number = self.root.get("number")
		self.title  = self.root.get("title") 
		self.link   = self.root.get("link")

