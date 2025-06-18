import torch
import torch.nn.functional as F
from torch import nn
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.data import Data


class QNetwork(nn.Module):
    """
    Graph Neural Network (GNN) based Q-Network for reinforcement learning. It is intended to work
    as a baseline and no bond features are used in the GNN layers.
    The network consists of two GCN layers followed by a multi-layer perceptron (MLP) head.
    """

    def __init__(self, in_features : int, hidden_features : int, num_actions : int):
        super().__init__()
        #GNN encoder
        self.conv1 = GCNConv(in_features, hidden_features)
        self.conv2 = GCNConv(hidden_features, hidden_features)

        #MLP head
        self.lin1 = nn.Linear(hidden_features, hidden_features)
        self.lin2 = nn.Linear(hidden_features, num_actions)

    
    def forward(self, graph_data : Data):

        x, edge_index = graph_data.x, graph_data.edge_index

        # GNN layers
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        x = F.relu(x)

        # Global pooling
        x = global_mean_pool(x, graph_data.batch)

        # MLP head
        x = self.lin1(x)
        x = F.relu(x)
        x = self.lin2(x)

        return x