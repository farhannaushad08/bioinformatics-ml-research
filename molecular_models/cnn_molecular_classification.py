"""
Classical CNN models for molecular data (prototype)

This script trains and evaluates on biomedical images . Replace with your biomedical images as needed.

"""

import tensorflow as tf
from tensorflow.keras import layers, models, datasets

(x_train, y_train), (x_test, y_test) = datasets.cifar10.load_data()

# Keep only 2 classes (binary example)
x_train, y_train = x_train[y_train.flatten() < 2], y_train[y_train.flatten() < 2]
x_test, y_test = x_test[y_test.flatten() < 2], y_test[y_test.flatten() < 2]

x_train, x_test = x_train / 255.0, x_test / 255.0

model = models.Sequential([
    layers.Conv2D(32, (3,3), activation='relu', input_shape=(32,32,3)),
    layers.MaxPooling2D((2,2)),
    layers.Conv2D(64, (3,3), activation='relu'),
    layers.MaxPooling2D((2,2)),
    layers.Flatten(),
    layers.Dense(64, activation='relu'),
    layers.Dense(2, activation='softmax')
])

model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
history = model.fit(x_train, y_train, epochs=5, validation_data=(x_test, y_test))
