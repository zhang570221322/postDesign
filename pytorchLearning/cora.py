import os
import torch
import torch.nn as nn
import torch.nn.functional as F
# PyG处理好的一些数据，如"Cora", "CiteSeer" and "PubMed" ，用Planetoid这个类调用即可
from torch_geometric.datasets import Planetoid
import torch_geometric.nn as pyg_nn


# 第一步：准备数据，load Cora dataset
def get_data(folder="node_classify/cora", data_name="cora"):
    """
    :param folder:保存数据集的根目录。
    :param data_name:数据集的名称
    :return:返回的是一个对象，就是PyG文档里的Data对象，它有一些属性，如 data.x、data.edge_index等
    """
    dataset = Planetoid(root=folder, name=data_name)
    return dataset

# 第二步：定义模型，create the graph cnn model


class GraphCNN(nn.Module):
    def __init__(self, in_c, hid_c, out_c):
        super(GraphCNN, self).__init__()  # 表示子类GraphCNN继承了父类nn.Module的所有属性和方法.
        # 下面这个就是前面讲的GCN，参数只有输入和输出，定义了两层的GCN.
        self.conv1 = pyg_nn.GCNConv(in_channels=in_c, out_channels=hid_c)
        self.conv2 = pyg_nn.GCNConv(in_channels=hid_c, out_channels=out_c)

    def forward(self, data):
        # data.x  data.edge_index
        x = data.x  # [N, C], C为特征的维度
        edge_index = data.edge_index  # [2, E], E为边的数量
        # [N, D], N是节点数量，D是第一层输出的隐藏层的维度
        hid = self.conv1(x=x, edge_index=edge_index)
        hid = F.relu(hid)
        # [N, out_c], out_c就是定义的输出，比如分成几类就是几，这里是7
        out = self.conv2(x=hid, edge_index=edge_index)
        out = F.log_softmax(out, dim=1)  # [N, out_c],表示输出
        return out


# todo list
class YouOwnGCN(nn.Module):  # 这个不用理会，如果之后想用别的图卷积实现，可以自己在这里写，然后调用
    pass


os.environ["CUDA_VISIBLE_DEVICES"] = "0"  # 配置GPU
cora_dataset = get_data()
# todo list
# 这个是自己写的网络的实例化
my_net = GraphCNN(in_c=cora_dataset.num_node_features,
                  hid_c=13, out_c=cora_dataset.num_classes)
device = torch.device(
    "cuda" if torch.cuda.is_available() else "cpu")   # 检查设备
my_net = my_net.to(device)  # 模型送入设备
data = cora_dataset[0].to(device)  # 数据送入设备，也就是一张图
# 第三步：定义损失函数和优化器
optimizer = torch.optim.Adam(my_net.parameters(), lr=1e-3)  # 优化器
# 第四步：训练+测试
# model train,这个train就是说归一化等可以重复使用，而设置成eval则就不行了，表示测试
my_net.train()
for epoch in range(200):
    optimizer.zero_grad()  # 每次缓存之后清零,不然梯度会累加
    output = my_net(data)  # 预测结果
    loss = F.nll_loss(output[data.train_mask],
                      data.y[data.train_mask])  # 意思就是只取训练集
    loss.backward()
    print("epoch:", epoch + 1, loss.item())
    optimizer.step()  # 优化器
# model test
my_net.eval()
_, prediction = my_net(data).max(dim=1)
target = data.y
test_correct = prediction[data.test_mask].eq(
    target[data.test_mask]).sum().item()
test_number = data.test_mask.sum().item()
print("Accuracy of Test Samples:{}%".format(100*test_correct/test_number))
