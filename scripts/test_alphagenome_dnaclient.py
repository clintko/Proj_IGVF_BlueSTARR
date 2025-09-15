### set environment
import numpy  as np
import pandas as pd
import os

from alphagenome.data   import genome
from alphagenome.models import dna_client, variant_scorers

### check API key exist
assert "ALPHAGENOME_API_KEY" in os.environ

### load alphagenome model
dna_model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])
print(dna_model)