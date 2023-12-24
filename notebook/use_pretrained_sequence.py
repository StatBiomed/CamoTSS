import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from hilearn.plot import ROC_plot
from sklearn.metrics import confusion_matrix

os.environ['CUDA_VISIBLE_DEVICES'] = '0'
## configure other parameter，如batch_size, num_workers, learning rate, 以及总的epochs
batch_size = 256
num_workers = 4   
lr = 1e-4
epochs = 2

allfeaturedf=pd.read_csv('/data/users/houruiyan/sequence_CNN/four_sequence_feature.csv',index_col=0)
allfeaturedf

train_df, test_df = train_test_split(allfeaturedf, test_size=0.2,random_state=666,shuffle=True)
train_df,validation_df=train_test_split(train_df,test_size=0.25,random_state=666,shuffle=True)

len(test_df)


class FMDataset(Dataset):
    def __init__(self, df, transform=None):
        self.df = df
        self.transform = transform
        self.images = df.iloc[:,7:].values.astype(np.uint8)
        self.labels = df['label'].values
        self.extrafeatures=df.iloc[:,3:7].values
        
    def __len__(self):
        return len(self.images)
    
    def __getitem__(self, idx):
        image = self.images[idx].reshape(4,200)
        #print(idx)
        
        label = int(self.labels[idx])
        extrafeature=self.extrafeatures[idx]
        
        if self.transform is not None:
            image = self.transform(image)
        else:
            image = torch.tensor(image, dtype=torch.float)
            
        label = torch.tensor(label, dtype=torch.long)
        extrafeature=torch.tensor(extrafeature,dtype=torch.float)
        
        return image, label,extrafeature

 

train_data = FMDataset(train_df)
test_data = FMDataset(test_df)
validation_data=FMDataset(validation_df)

train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True, num_workers=num_workers,drop_last=True)
test_loader = DataLoader(test_data, batch_size=batch_size, shuffle=False, num_workers=num_workers)
validation_loader=DataLoader(validation_data, batch_size=batch_size, shuffle=False, num_workers=num_workers)

len(test_loader)

import matplotlib.pyplot as plt
image, label,extrafeature = next(iter(test_loader))
print(image.shape, label.shape,extrafeature.shape)
plt.imshow(image[0], cmap="gray")

## load model

class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(4, 128, 8),
             nn.ReLU(),
            nn.Conv1d(128,64,4),
            nn.ReLU(),
             nn.BatchNorm1d(64),
              nn.MaxPool2d(2),
            nn.Dropout(0.4),
             nn.Flatten(),
            nn.Linear(3040, 32),
            nn.ReLU(),
        )
        self.fc = nn.Sequential(
            nn.Linear(32, 2),
            nn.Sigmoid()
        )
        
    def forward(self, x):
        x = self.conv(x)
        #print(x.shape)
        #combinex=torch.cat((x,extrafeature),dim=1)
        #x = x.view(-1, 3040)
        x = self.fc(x)
        # x = nn.functional.normalize(x)
        return x
    
    

    

model = Net()
model = model.cuda()

check_point=torch.load('/data/users/houruiyan/sequence_CNN/best_model.pth')
check_point

model.load_state_dict(check_point['model_state_dict'])

criterion = check_point['loss']

## Do test

def test(epoch):       
    model.eval()
    test_loss = 0
    gt_labels = []
    pred_labels = []
    y_scorels=[]
    with torch.no_grad():
        for data, label,extrafeature in test_loader:
            #print(data)
            data, label,extrafeature = data.cuda(), label.cuda(),extrafeature.cuda()
            output = model(data)
            #print(output.shape)
            preds = torch.argmax(output, 1)
            gt_labels.append(label.cpu().data.numpy())
            pred_labels.append(preds.cpu().data.numpy())
            loss = criterion(output, label)
            test_loss += loss.item()*data.size(0)
            y_scorels.append(output[:,1].cpu().data.numpy())
            
            #y_scorels.extend(y_score.tolist())
            
    # test_loss = test_loss/len(test_loader.dataset)
    # gt_labels, pred_labels,y_scorels = np.concatenate(gt_labels), np.concatenate(pred_labels),np.concatenate(y_scorels)
    # acc = np.sum(gt_labels==pred_labels)/len(pred_labels)
    # auc=roc_auc_score(gt_labels, y_scorels)
    # confu=confusion_matrix(gt_labels, pred_labels)
    
    
    #print('Epoch: {} \tValidation Loss: {:.6f}, Accuracy: {:6f},AUC:{:6f}'.format(epoch, val_loss, acc,auc))
    #print(confu)
    return gt_labels,y_scorels

gt_labels,y_scorels=test(1)

len(gt_labels)

20*256

groudtruthlabel=np.concatenate(gt_labels)
groudtruthlabel

y_scorels=np.concatenate(y_scorels)
y_scorels

res = ROC_plot(groudtruthlabel, y_scorels, threshold=0.5)

gtls=list(groudtruthlabel)
gtls

yscolels=list(y_scorels)
yscolels

onlysequencedf=pd.DataFrame({'Oly_sequence_groudTruth':gtls,'Oly_sequence_Score':yscolels})
onlysequencedf

onlysequencedf.to_csv('/data/users/houruiyan/sequence_CNN/only_sequence_Rscore.csv')
