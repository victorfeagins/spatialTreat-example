from datetime import datetime, timedelta

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from math import sin, cos, pi  # rotating regions
from math import floor  # truncating naics codes
import numba  # speed up data transform with JIT compilation

print("Did this run?")
