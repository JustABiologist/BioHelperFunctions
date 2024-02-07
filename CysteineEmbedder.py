import numpy as np
import torch
from torch import nn, optim
from torch.utils.data import DataLoader, TensorDataset, Dataset
from sklearn.model_selection import train_test_split
import pickle

def load_dict_from_pkl(filepath):
    """
    Loads and returns a dictionary from a pickle file.

    Args:
    - filepath (str): The path to the .pkl file.

    Returns:
    - dict: The dictionary loaded from the pickle file.
    """
    with open(filepath, "rb") as file:
        data_dict = pickle.load(file)
    return data_dict

def normalize_and_pad_data(data_dict, target_length=256):
    """
    Normalizes the data in the dictionary, pads arrays to a specific length, and flattens them.

    Args:
    - data_dict (dict): The dictionary containing the data to be normalized and padded.
    - target_length (int): The target length to pad the arrays to.

    Returns:
    - dict: A new dictionary with normalized, padded, and flattened data.
    """
    normalized_padded_dict = {}
    max_value = max(max(np.array(value).flatten()) for value in data_dict.values())
    
    for key, value in data_dict.items():
        # Normalize
        normalized_array = np.array(value) / max_value
        
        # Flatten
        flattened_array = normalized_array.flatten()
        if len(flattened_array) > target_length:
            continue
        
        # Pad
        pad_length = max(0, target_length - len(flattened_array))
        padded_array = np.pad(flattened_array, (0, pad_length), 'constant', constant_values=0)
        
        normalized_padded_dict[key] = (padded_array, padded_array)
    
    return normalized_padded_dict

def split_data(data_dict, test_size=0.1, val_size=0.15):
    """
    Splits the dictionary into training, validation, and test sets based on the specified sizes.
    """
    # Split keys for train and test
    keys = list(data_dict.keys())
    train_keys, test_keys = train_test_split(keys, test_size=test_size, random_state=42)
    
    # Further split train keys for validation
    train_keys, val_keys = train_test_split(train_keys, test_size=val_size, random_state=42)
    
    train_dict = {key: data_dict[key] for key in train_keys}
    val_dict = {key: data_dict[key] for key in val_keys}
    test_dict = {key: data_dict[key] for key in test_keys}
    
    return train_dict, val_dict, test_dict

class CustomDataset(Dataset):
    """
    A custom dataset to handle data in the form of dictionaries.
    """
    def __init__(self, data_dict):
        self.data = list(data_dict.values())
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        input, target = self.data[idx]
        return input, target


class Autoencoder(nn.Module):
    def __init__(self, input_size):
        super(Autoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_size, input_size//2),
            nn.ReLU(),
            nn.Linear(input_size//2, input_size//4),
            nn.ReLU(),
            nn.Linear(input_size//4, 2)
        )
        self.decoder = nn.Sequential(
            nn.Linear(2, input_size//4),
            nn.ReLU(),
            nn.Linear(input_size//4, input_size//2),
            nn.ReLU(),
            nn.Linear(input_size//2, input_size),
            nn.ReLU()  # Assuming you want the ReLU activation here as well
        )

    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        return x


def train_model(model, train_dict, val_dict, test_dict, epochs=1000, batch_size=32, learning_rate=1e-3):
    """
    Trains the model using the provided training and validation dictionaries, and evaluates on the test dictionary.
    """
    # Convert dictionary values to tensors within the CustomDataset
    print(len(train_dict))
    train_loader = DataLoader(CustomDataset(train_dict), batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(CustomDataset(val_dict), batch_size=batch_size, shuffle=False)
    test_loader = DataLoader(CustomDataset(test_dict), batch_size=batch_size, shuffle=False)
    
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        for inputs, targets in train_loader:
            inputs, targets = torch.tensor(inputs, dtype=torch.float32), torch.tensor(targets, dtype=torch.float32)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * inputs.size(0)
        
        train_loss /= len(train_loader.dataset)

        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for inputs, targets in val_loader:
                inputs, targets = torch.tensor(inputs, dtype=torch.float32), torch.tensor(targets, dtype=torch.float32)
                outputs = model(inputs)
                loss = criterion(outputs, targets)
                val_loss += loss.item() * inputs.size(0)
        
        val_loss /= len(val_loader.dataset)

        print(f'Epoch {epoch+1}: Train Loss: {train_loss:.4f}, Validation Loss: {val_loss:.4f}')

    # Evaluate on test set and return outputs in a dictionary
    test_outputs = {}
    model.eval()
    with torch.no_grad():
        for idx, (inputs, targets) in enumerate(test_loader):
            inputs = torch.tensor(inputs, dtype=torch.float32)
            outputs = model(inputs)
            key = str(idx)  # Simple unique identifier for each input
            test_outputs[key] = outputs.numpy().tolist()
    
    return model, test_outputs

def main(input_path):
    # Load the data dictionary from a pickle file
    data_dict = load_dict_from_pkl(input_path)

    # Normalize and pad the data
    normalized_padded_dict = normalize_and_pad_data(data_dict)

    # Split the data into training, validation, and test sets
    train_dict, val_dict, test_dict = split_data(normalized_padded_dict)

    # Prepare the input size for the autoencoder model
    # Assuming all data points in normalized_padded_dict have the same length after padding
    input_size = len(next(iter(normalized_padded_dict.values()))[0])

    # Initialize the autoencoder model
    model = Autoencoder(input_size)

    # Train the model
    model, test_outputs = train_model(model, train_dict, val_dict, test_dict)

    # Here you can save the trained model and test_outputs if needed
    print("Training completed.")


input_path = '/Users/floriangrun/Desktop/data_cysteine_dmaps.pkl'
main(input_path)