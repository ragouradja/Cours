library(keras)
install_keras()
train=sample(1:150,100)
x_train = as.matrix(iris[train, 1:4])
y_train = iris[train, 5]
x_test = as.matrix(iris[-train, 1:4])
y_test = iris[-train,5]
y_train = to_categorical(as.integer(y_train)-1)
y_test = to_categorical(as.integer(y_test)-1)
model = keras_model_sequential()
layer_dense(model, units=5, activation="sigmoid", input_shape=4)
layer_dense(model, units=3, activation="softmax")
compile(model, loss = "categorical_crossentropy", optimizer = optimizer_adam(), metrics = "accuracy")
history = fit(model, x_train, y_train, epochs = 100, batch_size = 20, validation_split = 0.2)

classes=predict_classes(model, x_train)
table(classes, iris[train,5])
classes=predict_classes(model, x_test)
table(classes, iris[-train,5])
evaluate(model, x_test, y_test)

