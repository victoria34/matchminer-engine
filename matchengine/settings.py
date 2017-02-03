"""Copyright 2016 Dana-Farber Cancer Institute"""

import os
import sys
import json
import logging

# logging
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(message)s', )

TUMOR_TREE = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data/tumor_tree.txt'))

months = [
    'January', 'February', 'March', 'April', 'May', 'June',
    'July', 'August', 'September', 'October', 'November', 'December'
]

MONGO_URI = "mongodb://localhost:27017/matchminer?replicaSet=rs0"

uri_check = os.getenv("MONGO_URI", None)
if uri_check:
    MONGO_URI = uri_check