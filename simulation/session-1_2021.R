# Description: Session 1 live
# Date: 2021-09-20
# Author: BT


require(MASS)  # pour rendre disponible la fonction "truehist"



#-----------
# Aiguille de Buffon

rm(list = ls())  # pour partir sur des bases saines: effacer la mémoire


# Estimateur de pi par Aiguille de Buffon
# r: rapport longueur aiguille / espace inter-lame de parquet
# n: nombre d'aiguilles lancées
estim_pi <- function(n, r) {
	theta <- runif(n, min = 0, max = pi) # angle aléatoire: unif sur [0, pi]
	x <- runif(n, min = 0, max = 1)  # position de la tete: unif entre [0 et 1]

	a <- x + cos(theta) * r  # position de la pointe de l'aiguille

	ns <- sum(a > 1 | a < 0)  # nombre d'aiguilles entre deux lames
	
	ps <- ns / n  # proportion...
	
	return(2 * r / ps)
}


estim_pi(n = 100, r = 0.3)  # Une seule estimation de pi



# Accumuler nrep estimations de pi
mse_pi <- function(n, r, nrep) {
	
	pi_hat <- rep(NA, times = nrep)  # préparation d'un espace suffisant
	for (i in 1:nrep) {  # accumulation
		pi_hat[i] <- estim_pi(n, r)
	}

	mse <- mean((pi_hat - pi)^2)

	
	truehist(pi_hat, xlim = c(0, 8))  # histogramme
	abline(v = pi, col = "red")  # vraie valeur
	title(main = paste("mse = ", round(mse, 6)))
	
	
	return(mse)
}


mse_pi(n = 100, r = 0.3, nrep = 100000)  # caractériser l'estimaeur (n = 100, r = 0.3) à l'aide de 100000 répétitions



#-----------
# Fonction d'une variable aléatoire

rm(list = ls())  # pour repartir sur des bases saines: effacer la mémoire


nrep <- 100000  # pour bien caractériser la variable aléatoire produite, on fait 100000 tirages

# echantillonnage dans une loi normale centré-reduite "N(0, 1)"
x <- rnorm(n = nrep, mean = 0, sd = 1)

require(MASS)
truehist(x)  # pour vérifier que la densité ressemble courbe en cloche

# surimposition de la densité d'une loi normale CR
val_x <- seq(from = -4, to = 4, len = 1000)  # NB: j'ai augmenté le nombre de points pour construire la courbe
val_y <- dnorm(val_x)
lines(val_x, val_y, col = "red", lwd = 2)  
# On s'attend à ce que, si la fonction "rnorm" a étyé correctement implémentée, à ce que
# la courbe rouge de la fonction de densité "théorique" suive très exactement les contours de l'histogramme




# estimation de l'espérance
mean(x)  # on doit avoir une valeur proche de 0 car l'espérance d'une loi N(0, 1) est 0

# estimation de la variance
var(x)  # on doit avoir une valeur proche de 1 car l'espérance d'une loi N(0, 1) est 1


# estimation du quantile 2.5% et 97.5%
quantile(x, probs = c(0.025, 0.975))  # on doit avoir des valeurs proches de -1.959964 et +1.959964



# Exemple de fonction de variable aléatoire: carré
# Variable: Y = X²
y <- x^2

# Distribution de X² (avec X normal CR)
truehist(y, xlim = c(0, 7))
# La courbe est en "L"


# surimposition de la densité d'une loi de chi2 à 1 ddl (Pearson a en effet donné le nom de "chi2" à la loi du carré d'une loi normale centrée-réduite)
val_x <- seq(from = 0, to = 7, len = 1000)
val_y <- dchisq(val_x, df = 1)
lines(val_x, val_y, col = "red", lwd = 2)





# On fait tente maintenant le carré d'une loi uniforme
# echantillonnage dans une loi uniforme
x <- runif(n = nrep)

# histogramme pour visualiser la distribution
truehist(x)  # pour vérifier qu'on part bien d'une loi "plate"

# On prend maintenant le carré
# Y = X² (mais cette fois-ci X suit une loi uniforme et pas une loi normale)
y <- x^2

# Distriubution de X²
truehist(y)
# c'est une loi en "L" 

# surimposition de la densité d'une loi de chi2 à 1 ddl
val_x <- seq(from = 0, to = 7, len = 1000)
val_y <- dchisq(val_x, df = 1)
lines(val_x, val_y, col = "red", lwd = 2)

# mais qui ne correspond plus à une loi de chi² à 1 ddl




# Autre transformation: transformation linéaire
# Y = X + b X


# echantillonnage dans une loi uniforme
x <- runif(n = nrep, min = 0, max = 1)

range(y) # support [0, 1]
round(range(y), 5)

# distribution: plate
truehist(x)
val_x <- seq(from = 0, to = 1, len = 100)
val_y <- dunif(val_x, min = 0, max = 1)
lines(val_x, val_y, col = "red", lwd = 2)


# On choisit: a = -3 (translation) et b = 2 (mise à l'échelle)
y <- -3 + 2*x

range(y) # on a maintenant un support [-3, -1] = a + b *[0, 1] = -3 * 2 * [0, 1] = -3 + [0, 2]
round(range(y), 5)

# distribution: toujours "plate" -> toujours uniforme
truehist(y)

val_x <- seq(from = -3, to = -1, len = 100)
val_y <- dunif(val_x, min = -3, max = -1)  # mais celle d'une loi uniforme sur [-3, -1]
lines(val_x, val_y, col = "red", lwd = 2)


# Alternative: slide 38
a <- -3
b <- 2
val_y <- dunif((val_x - a) / b, min = 0, max = 1) / abs(b)
lines(val_x, val_y, col = "blue", lwd = 2)




#-----------
# INTEGRATION

rm(list = ls())  # pour repartir sur des bases saines: effacer la mémoire


# fonction à intégrer
g <- function(x) {
	return((exp(x) - 1) / (exp(1) - 1))		
}
	


# graphe de g
val_x <- seq(from = 0, to = 2, len = 100)
val_y <- g(val_x)
plot(val_x, val_y, type = "l")


# majorant choisi: g(2)
abline(h = g(2), col = "red")  # on peut le visualiser sur le graphique


# I attendu: calcul analytique
I_expected <- (exp(2) - 3)/(exp(1) - 1)
print(I_expected)

# On peut vérifier à l'aide de la fonction 'integrate'
integrate(g, lower = 0, upper = 2)



# Estimateur de I par la methode B&W
estim_I_BW <- function(n) {
	b <- 2
	a <- 0
	
	# majorant de g: on peut prendre g(2)
	# mais on aurait aussi pu prendre une valeur plus grande, par exemple 5, mais ça produirait un autre estimateur avec des propriétés différentes
	m <- g(2)

	u <- runif(n, min = a, max = b)  # abscisse aléatoire
	v <- runif(n, min = 0, max = m)  # ordonnée aléatoire

	ns  <- sum(v < g(u))  # nb de points aleatoires sous la courbe "g"
	
	I_estim <- m * (b - a) * ns / n 
		
	return(I_estim)
}


estim_I_BW(100)  # Une seule estimation de I par la méthode B&W avec n = 100



# Accumuler nrep estimations de I par la méthode B&W
mse_I_BW <- function(n, nrep) {
	
	I_hat <- rep(NA, times = nrep)
	for (i in 1:nrep) {
		I_hat[i] <- estim_I_BW(n)
	}

	mse <- mean((I_hat - I_expected)^2)

	
	truehist(I_hat, xlim = c(1, 5))  # histogramme
	abline(v = I_expected, col = "red")  # vraie valeur
	title(main = paste("BW - n = ", n, "\nmse = ", round(mse, 6)))
	
	
	return(mse)
}

mse_I_BW(n = 100, nrep = 100000)  # pour caractériser cet estimateur B&W avec n = 100




# Echantillonnage simple

# Estimateur de I par la methode d'échantillonnage Simple
estim_I_simple <- function(n) {

	b <- 2
	a <- 0
	
	# echantillonnage uniforme entre a et b
	x <- runif(n, min = a, max = b)
	
	# "I = (b-a) E[g(X)]"
	# "I_estime = (b-a) * sum(g(xi))/n = (b-a) moyenne(g(xi))"
	I_estim <- (b - a) * mean(g(x))
	
	return(I_estim)
}


estim_I_simple(100)  # Une seule estimation de I par la méthode d'échantillonnage simple avec n = 100



# accumulons...
mse_I_simple <- function(n, nrep) {
	
	I_hat <- rep(NA, times = nrep)
	for (i in 1:nrep) {
		I_hat[i] <- estim_I_simple(n)
	}

	mse <- mean((I_hat - I_expected)^2)

	
	truehist(I_hat, xlim = c(1, 5))  # histogramme
	abline(v = I_expected, col = "red")  # vraie valeur
	title(main = paste("Simple - n = ", n, "\nmse = ", round(mse, 6)))
	
	
	return(mse)
}

mse_I_simple(n = 100, nrep = 100000)  # caractérisons

