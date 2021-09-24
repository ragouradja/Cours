library(keras)
install_keras()
mnist <- dataset_mnist()
rotate <- function(x) t(apply(x, 2, rev))
image(1:28,1:28,rotate(mnist$train$x[1,,]),col=grey.colors(12))
mnist$train$y[1]
x_train = mnist$train$x
y_train = mnist$train$y
x_test = mnist$test$x
y_test = mnist$test$y

# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# rescale
x_train <- x_train / 255
x_test <- x_test / 255
# output generation
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)


model = keras_model_sequential()
#784 car à plat
layer_dense(model, units=128, activation="relu", input_shape=784)
layer_dropout(model, rate = 0.4) # 40% à 0
layer_dense(model, units=64, activation="relu") # pas d'input shape car deja une couche avant
layer_dropout(model, rate = 0.3)

layer_dense(model, units=10, activation="softmax")

model = compile(model, loss = "categorical_crossentropy", optimizer = optimizer_adam(), metrics = "accuracy")
history = fit(model, x_train, y_train, epochs = 30, batch_size = 10, validation_split = 0.2)

classes=predict_classes(model, x_train)
table(classes, iris[train,5])
classes=predict_classes(model, x_test)
table(classes, iris[-train,5])
evaluate(model, x_test, y_test)

