p_0 <- 0.2

log_0_prepering <- function(p = 0.5){
  tdistr <- (1 - p_0) * (-p) / log(1 - p)
  sum_distr <- p_0
  sum <- c(p_0)
  k <- 2
  
  while(tdistr > 1e-10){
    sum_distr <- sum_distr + tdistr
    sum <- c(sum, sum_distr)
    tdistr <- tdistr * p / k * (k - 1)
    k <- k + 1
  }
  
  return(sum)
}

rlog_0 <- function(p = 0.5, sum){
  x <- runif(1)
  j <- 1
  
  while (x > sum[j]){
    j <- j + 1
  }
  
  return(j - 1)
}

rnlog_0 <- function(n = 1, p = 0.5, sum){
  v <- replicate(n, rlog_0(p, sum))
  
  return(v)
}

table_log <- log_0_prepering(0.75)
v <- rnlog_0(10000, 0.75, table_log)
hist(v[v < 9], probability = TRUE, breaks = seq(-0.5, max(v) + 0.5, 1), xlim = c(-0.5, 8.5), xlab = "Sample", main = "")

# Табличная функция распредления логарифмического распределения
log_prepering <- function(p = 0.5){
  tdistr <- (-p) / log(1 - p)
  sum_distr <- tdistr
  sum <- c(sum_distr)
  k <- 2
  
  while(tdistr > 1e-10){
    tdistr <- tdistr * p / k * (k - 1)
    sum_distr <- sum_distr + tdistr
    sum <- c(sum, sum_distr)
    k <- k + 1
  }
  
  return(sum)
}

# Моделирование логарифмического распределения по табличной функции распределения
rlog <- function(sum){
  x <- runif(1)
  j <- 1
  
  while (x > sum[j]){
    j <- j + 1
  }
  
  return(j)
}

# Моделирование выборки логарифмического распределения
rnlog <- function(n = 1, sum){
  v <- replicate(n, rlog(sum))
  
  if(n == 0)
    v <- 0;
  
  return(v)
}

# Моделирование биномиально-логарифмического распределения.
# sumlog - накопленные вероятности логарифмического распределения.
rbinomlog <- function(num = 1, n = 1, p = 0.5, sumlog){
  res <- c()
  
  for(i in rbinom(num, n, p)){
    res <- c(res, sum(0, rnlog(i, sumlog)))
  }
  
  return(res)
}

# Получение чисел Стирлинга первого рода
MakeNumStir <- function(data, n){
  for (i in 2:n){
    for (j in 2:n){
      if (i == j){
        data[i, j] <- 1 / gamma(i)
      } else {
        data[i, j] <- data[i - 1, j - 1] / (i - 1) + data[i - 1, j] * (i - 2) / (i - 1) 
      }
    }
  }
  
  return(data)
}

nNumStir <- 600
NumStir <- as.data.frame(matrix(nrow = nNumStir, ncol = nNumStir))
for (i in 1:nNumStir){
  NumStir[1, i] <- 0
  NumStir[i, 1] <- 0
}
NumStir[1, 1] <- 1
NumStir <- MakeNumStir(NumStir, nNumStir)

# Вероятности биномиально-логарифмического распределения
get_prob_binomlog <- function(k, p = 0.5, n = 1, q = 0.5){
  res <- 0
  alpha <- -1 / log(1 - q)
  frac <- 1

  for (j in 0:k){
    res <- res + frac * NumStir[k + 1, j + 1] * (p * alpha) ^j * (1 - p) ^(k - j)
    frac <- frac * (n - j)
  }
  
  prob <- -(1 - p) ^(n - k) * q ^k * res
  
  return((1 - p) ^(n - k) * q ^k * res)
}

# Функция правдоподобия для бином-логарифма
m_log_lik_binomlog <- function(x.in, p = 0.5, n = 1, q = 0.5){
  if (q <= 0 || q >= 1 || p <= 0 || p >= 1 || n <= 0){
    return(-log(0));
  }
  res <- 0

  for (i in x.in){
    prob <- get_prob_binomlog(i, p, n, q)
    
    if (prob < 0 || prob > 1 || is.nan(prob))
      return(-log(0))
    
    res <- res + log(prob)
  }

  return(-res)
}

# Моделирование логарифмически-биномиального распределения
rlogbinom <- function(num = 1, sumlog, n = 1, p = 0.5){
  res <- c()
  
  for(i in rnlog(num, sumlog)){
    res <- c(res, sum(0, rbinom(i, n, p)))
  }
  
  return(res)
}

# Вероятность лог-бинома для k = 0
get_prob_logbinom_logbinom_0 <- function(q = 0.5, p = 0.5, n = 1){
  return(log(1 - q* (1 - p) ^n) / log(1 - q))
}

# Вероятность лог-бинома для k = 1
get_prob_logbinom_logbinom_1 <- function(q = 0.5, p = 0.5, n = 1){
  return(-q * p * n * (1 - p)^(n - 1) / ((1 - q* (1 - p) ^n) * log(1 - q)))
}

# Вероятность бином-логарифма для k = 2
get_prob_logbinom_logbinom_2 <- function(q = 0.5, p = 0.5, n = 1){
  return(1 / 2 * (-q * p ^2 * n * (1 - p)^(n - 2) / ((1 - q* (1 - p) ^n) * log(1 - q))) *
           (q * n * (1 - p)^(n) / ((1 - q* (1 - p) ^n)) + (n - 1)))
}

# Вероятность лог-бинома для k = 3
get_prob_logbinom_logbinom_3 <- function(q = 0.5, p = 0.5, n = 1){
  return (1 / 6 * (-q * p ^3 * n * (1 - p)^(n - 3) / ((1 - q* (1 - p) ^n) * log(1 - q))) *
            (2 * q ^2 * n ^2 * (1 - p)^(2 * n) / ((1 - q* (1 - p) ^n) ^2) +
               3 * q * n * (n - 1) * (1 - p)^(n) / ((1 - q* (1 - p) ^n)) + (n - 1) * (n - 2)))
}

# Вероятность лог-бинома для k = 4
get_prob_logbinom_logbinom_4 <- function(q = 0.5, p = 0.5, n = 1){
  return (1 / 24 * (-q * p ^4 * n * (1 - p)^(n - 4) / ((1 - q* (1 - p) ^n) * log(1 - q))) *
            (3 * q ^3 * n ^3 * (1 - p)^(3 * n) / ((1 - q* (1 - p) ^n) ^3) +
               12 * q ^2 * n ^2 * (n - 1) * (1 - p)^(2 * n) / ((1 - q* (1 - p) ^n) ^2) +
               q * (1 - p)^(n) / ((1 - q* (1 - p) ^n)) * (3 * n * (n - 1) ^2 + 4 * n * (n - 1) * (n - 2)) +
               (n - 1) * (n - 2) * (n - 3)))
}

# Вероятности логарифмически-биномиального распределения
get_prob_logbinom <- function(k, q = 0.5, p = 0.5, n = 1){
  n <- ceiling(n)
  
  if (k == 0){
    return(get_prob_logbinom_logbinom_0(q, p, n))
  }
  
  if (k == 1){
    return(get_prob_logbinom_logbinom_1(q, p, n))
  }
  
  if (k == 2){
    return(get_prob_logbinom_logbinom_2(q, p, n))
  }
  
  if (k == 3){
    return(get_prob_logbinom_logbinom_3(q, p, n))
  }
  
  if (k == 4){
    return(get_prob_logbinom_logbinom_4(q, p, n))
  }
  
  if (k >= 5){
    return(1 - get_prob_logbinom_logbinom_0(q, p, n) - get_prob_logbinom_logbinom_1(q, p, n) -
           get_prob_logbinom_logbinom_2(q, p, n) - get_prob_logbinom_logbinom_3(q, p, n) -
           get_prob_logbinom_logbinom_4(q, p, n))
  }
}

# Функция правдоподобия для лог-бинома
m_log_lik_logbinom <- function(x.in, q = 0.5, p = 0.5, n = 1){
  if (q <= 0 || q >= 1 || p <= 0 || p >= 1 || n <= 0){
    return(-log(0));
  }
  res <- 0
  
  for (i in x.in){
    prob <- get_prob_logbinom(i, q, p, n)

    if (prob < 0 || prob > 1 || is.nan(prob))
      return(-log(0))

    res <- res + log(prob)
  }
  
  return(-res)
}

# Моделирование свёртки отрицательно-биномиальных распределений
rdoublenbinom <- function(num = 1, p_1 = 0.5, size_1= 1, p_2 = 0.5, size_2 = 1){
  return(rnbinom(num, size_1, p_1) + rnbinom(num, size_2, p_2))
}

# Вероятности свёртки негативных биномов
get_prob_doublenbinom <- function(k, p_1 = 0.5, size_1= 1, p_2 = 0.5, size_2 = 1){
  res <- 0
  
  for (m in 0:k){
    res <- res + choose(m + size_1 - 1, m) * choose(k - m + size_2 - 1, k - m) * ((1 - p_1) ** m) * ((1 - p_2) ** (k - m))
  }
  
  return(res * (p_1 ** size_1) * (p_2 ** size_2))
}

# Функция правдоподобия свёртки негативных биномов
m_log_lik_doublenbinom <- function(x.in, p_1 = 0.5, size_1= 1, p_2 = 0.5, size_2 = 1){
  if (p_1 < 0 | p_1 > 1 | p_2 < 0 | p_2 > 1 | size_1 <= 0 | size_2 <= 0)
    return(-log(0))
  
  res <- 0
  
  for (i in x.in){
    prob <- get_prob_doublenbinom(i, p_1, size_1, p_2, size_2)
    
    if (prob < 0 || prob > 1 || is.nan(prob))
      return(-log(0))
    
    res <- res + log(prob)
  }
  
  -res
}

# Функция правдоподобия негативного бинома
m_log_lik_nbinom <- function(x.in, q = 0.5, n = 1){
  if (q <= 0 || q >= 1 || n <= 0){
    return(-log(0));
  }
  res <- 0
  
  for (i in x.in){
    prob <- dnbinom(i, n, q)
    
    if (prob < 0 || prob > 1 || is.nan(prob))
      return(-log(0))
    
    res <- res + log(prob)
  }
  
  return(-res)
}

# Создание выборки по частотам значений
generate_sample <- function(num_k){
  sam <- c()
  
  for (i in 1:length(num_k)){
    sam <- c(sam, rep.int(i - 1, num_k[i]))
  }
  
  return(sam)
}

# Критерий хи-квадрат
my_chisq <- function(exp_prob, prob){
  res <- 0
  
  for (i in 1:length(prob)){
    res <- res + (exp_prob[i] - prob[i])^2/prob[i]
  }
  
  return(res)
}

# Проверка гипотезы о соответствии теоритического распредления эмперическому
hist_make <- function (n, exp_prob, get_hist, distr = 'binomlog'){
  sam <- generate_sample(exp_prob)
  N <- length(sam)
  var_sam = var(sam) * (N - 1) / N
  mean_sam = mean(sam)
  exp_prob <- exp_prob / N
  
  if (n == 0){
    df <- -1
  }
  else{
    df <- 0
  }

  if (distr == 'binomlog'){
    if (n == 0){
      left <- 1
      right <- 1000
      while (right - left > 2) {
        nLeft <- (2 * left + right) %/% 3
        nRight <- (left + 2 * right) %/% 3
        qLeft <- optimize(function(x) m_log_lik_binomlog(sam, p = (var_sam / mean_sam * log(1 - x) * (1 - x) - log(1 - x)) / x, n = nLeft, q = x), c(0.01, 0.99))
        qRight <- optimize(function(x) m_log_lik_binomlog(sam, p = (var_sam / mean_sam * log(1 - x) * (1 - x) - log(1 - x)) / x, n = nRight, q = x), c(0.01, 0.99))

        if (qLeft$objective <= qRight$objective) {
          right <- nRight
        }
        else{
          left <- nLeft
        }
      }
      
      minValue <- optimize(function(x) m_log_lik_binomlog(sam, p = (var_sam / mean_sam * log(1 - x) * (1 - x) - log(1 - x)) / x, n = left, q = x), c(0.01, 0.99))
      q <- minValue$minimum
      n <- left
      for (i in (left + 1):right){
        value <- optimize(function(x) m_log_lik_binomlog(sam, p = (var_sam / mean_sam * log(1 - x) * (1 - x) - log(1 - x)) / x, n = i, q = x), c(0.01, 0.99))
        if (value$objective <= minValue$objective){
          minValue <- value
          q <- value$minimum
          n <- i
        }
      }
    }
    else {
      q <- optimize(function(x) m_log_lik_binomlog(sam, p = (var_sam / mean_sam * log(1 - x) * (1 - x) - log(1 - x)) / x, n = n, q = x), c(0.01, 0.99))$minimum
    }

    res <- c(n, q, (var_sam / mean_sam * log(1 - q) * (1 - q) - log(1 - q)) / q)
  }
  else if (distr == 'logbinom'){
    if (n == 0){
      left <- 1
      right <- 1000
      while (right - left > 1) {
        nLeft <- (2 * left + right) %/% 3
        nRight <- (left + 2 * right) %/% 3
        qLeft <- optimize(function(x) m_log_lik_logbinom(sam, q = x, p = - log(1 - x) * (1 - x) / nLeft / x * mean_sam, n = nLeft), c(0.01, 0.99))
        qRight <- optimize(function(x) m_log_lik_logbinom(sam, q = x, p = - log(1 - x) * (1 - x) / nRight / x * mean_sam, n = nRight), c(0.01, 0.99))
        
        if (qLeft$objective <= qRight$objective) {
          right <- nRight - 1
        }
        else{
          left <- nLeft + 1
        }
      }
      
      minValue <- optimize(function(x) m_log_lik_logbinom(sam, q = x, p = - log(1 - x) * (1 - x) / left / x * mean_sam, n = left), c(0.01, 0.99))
      q <- minValue$minimum
      n <- left
      for (i in (left + 1):right){
        value <- optimize(function(x) m_log_lik_logbinom(sam, q = x, p = - log(1 - x) * (1 - x) / i / x * mean_sam, n = i), c(0.01, 0.99))
        if (value$objective <= minValue$objective){
          minValue <- value
          q <- value$minimum
          n <- i
        }
      }
    }
    else {
      q <- optimize(function(x) m_log_lik_logbinom(sam, q = x, p = - log(1 - x) * (1 - x) / n / x * mean_sam, n = n), c(0.01, 0.99))$minimum
    }
    
    res <- c(n, q, - log(1 - q) * (1 - q) / n / q * mean_sam)
  } else if (distr == 'nbinom'){
    q <- optim(c(0.2, 5), function(x) m_log_lik_nbinom(sam, x[1], x[2]))
    n <- q$par[2]
    q <- q$par[1]
    res <- c(n, q)
  } else if (distr == "doublenbinom"){
    # print((mean_sam / 0.5 - var_sam) * 0.8 / (1 - 0.8) * 0.8 * 0.5 * (0.8 - 0.5))
    # print((mean_sam / 0.8 - var_sam) * 0.5 / (1 - 0.5) * 0.5 * 0.8 * (0.5 - 0.8))
    q <- optim(c(0.8, 1, 0.5, 1),
               function(x) m_log_lik_doublenbinom(sam, p_1 = x[1],
                                                  size_1 = x[2],
                                                  p_2 = x[3],
                                                  size_2 = x[4]))
    res <- c(q$par[1], q$par[2], q$par[3], q$par[4])
    # res <- c(res[1], round((mean_sam / res[2] - var_sam) * res[1] / (1 - res[1]) * res[1] * res[2] * (res[1] - res[2])),
    #          res[2], round((mean_sam / res[1] - var_sam) * res[2] / (1 - res[2]) * res[2] * res[1] * (res[2] - res[1])))
  }
  
  # print(res)

  prob <- c()
  
  for (i in 0:(length(exp_prob) - 2)){
    if (distr == 'binomlog'){
      prob <- c(prob, get_prob_binomlog(i, p = res[3], n = res[1], q = res[2]))
    }
    else if (distr == 'logbinom'){
      prob <- c(prob, get_prob_logbinom(i, q = res[2], p = res[3], n = res[1]))
    }
    else if (distr == 'nbinom'){
      prob <- c(prob, dnbinom(i, size = res[1], prob = res[2]))
    }
    else if (distr == 'doublenbinom'){
      prob <- c(prob, get_prob_doublenbinom(i, p_1 = res[1], size_1 = res[2], p_2 = res[3], size_2 = res[4]))
    }
  }
  
  prob <- c(prob, 1 - sum(prob))
  
  max_el_x = max(sam)
  max_el_y = max(exp_prob, prob)
  
  if (get_hist){
    hist.default(sam, probability = TRUE, breaks = seq(-0.5, max_el_x + 0.5, 1), xlim = c(-0.5, max_el_x + 0.5), ylim = c(0, max_el_y), xlab = "Sample", main = "")
    points(0:(length(prob) - 1), prob)
  }
  
  for (i in 1:length(prob)){
    while (100 * prob[i] < 5 && i < length(prob)){
      prob[i] <- prob[i] + prob[i + 1]
      prob <- prob[-(i + 1)]
      exp_prob[i] <- exp_prob[i] + exp_prob[i + 1]
      exp_prob <- exp_prob[-(i + 1)]
    }
  }
  
  if (100 * prob[length(prob)] < 5){
    prob[length(prob) - 1] <- prob[length(prob) - 1] + prob[length(prob)]
    prob <- prob[-length(prob)]
    exp_prob[length(exp_prob) - 1] <- exp_prob[length(exp_prob) - 1] + exp_prob[length(exp_prob)]
    exp_prob <- exp_prob[-length(exp_prob)]
  }
  
  df <- df + length(prob) - 3

  if (df < 1){
    df <- 1
  }
  
  chi <- N * my_chisq(exp_prob, prob)
  res <- c(res, 1 - pchisq(chi, df))
  
  return(res)
}

# Функция правдоподобия для ЛБР и БЛР
funcMP_3D <- function(data, distr = 'binomlog', inter_par1 = c(0.01, 0.99, 100), inter_par2 = c(1, 100, 10)){
  library(rgl)
  
  open3d()
  lines3d(c(0, 0), c(0, 0), c(0, 12), color = "gray");
  lines3d(c(0, 12), c(0, 0), c(0, 0), color = "gray");
  lines3d(c(0, 0), c(0, 12), c(0, 0), color = "gray");
  
  sam <- generate_sample(data)
  N <- length(sam)
  var_sam = var(sam) * (N - 1) / N
  mean_sam = mean(sam)
  x <- c()
  y <- c()
  z <- c()
  for (par1 in seq(inter_par1[1], inter_par1[2], (inter_par1[2] - inter_par1[1]) / inter_par1[3])){
    for (par2 in seq(inter_par2[1], inter_par2[2], (inter_par2[2] - inter_par2[1]) / inter_par2[3])){
      if (distr == 'binomlog'){
        value <- m_log_lik_binomlog(sam, p = (var_sam / mean_sam * log(1 - par1) * (1 - par1) - log(1 - par1)) / par1, n = par2, q = par1)
      }
      else if (distr == 'logbinom'){
        value <- m_log_lik_logbinom(sam, q = par1, p = - log(1 - par1) * (1 - par1) / par2 / par1 * mean_sam, n = par2)
      }
      
      if (value != 0) {
        x <- c(x, par1)
        y <- c(y, par2)
        z <- c(z, value)
      }
    }
  }
  
  z <- z - min(z)
  
  if (distr == 'binomlog'){
    z <- sqrt(z) / 10
    points3d(10 * x[z < 12], y[z < 12] * 10 / inter_par2[2], z[z < 12], color ="black")
  }
  else if (distr == 'logbinom'){
    points3d(10 * x[z < 12], y[z < 12] * 10 / inter_par2[2], z[z < 12], color ="black")
  }
}

# Проверка оценки для ЛБР
samLBR <- rlogbinom(100, log_prepering(0.4), 6, 0.3)
samLBR <- sapply(samLBR, function(x) if (x > 4){5} else {x})
histLBR <- hist(as.numeric(samLBR), probability = TRUE, breaks = seq(-0.5, max(samLBR) + 0.5, 1), xlim = c(-0.5, max(samLBR) + 0.5), plot = TRUE)
res <- hist_make(0, as.numeric(histLBR$counts), get_hist = TRUE, 'logbinom')
funcMP_3D(histLBR$counts, distr = 'logbinom')
res

# Проверка оценки для БЛР
samBLR <- rbinomlog(100, 6, 0.3, log_prepering(0.4))
histBLR <- hist(as.numeric(samBLR), probability = TRUE, breaks = seq(-0.5, max(samBLR) + 0.5, 1), xlim = c(-0.5, max(samBLR) + 0.5), plot = TRUE)
res <- hist_make(0, as.numeric(histBLR$counts), get_hist = TRUE, 'binomlog')
funcMP_3D(histBLR$counts)
res

# Загрузка радиобиологических данных
m <- read.csv("./VitroVivo.csv")
df.BLR <- data.frame(n = rep(0, 19), q = rep(0, 19), p = rep(0, 19), p_value = rep(0, 19))
df.LBR <- data.frame(n = rep(0, 19), q = rep(0, 19), p = rep(0, 19), p_value = rep(0, 19))

# Нахождение n для которого p-value максимально
for (i in 1:19){
  df.BLR[i, ] <- hist_make(0, as.numeric(m[i,]), get_hist = FALSE, 'binomlog')
  df.LBR[i, ] <- hist_make(0, as.numeric(m[i,]), get_hist = FALSE, 'logbinom')

  print(i)
}

# Гистограммы для строки данных
res <- hist_make(0, as.numeric(m[14,]), get_hist = TRUE, 'binomlog')
res
res <- hist_make(0, as.numeric(m[14,]), get_hist = TRUE, 'logbinom')
res

# Функция правдоподобия
funcMP_3D(as.numeric(m[1,]), distr = 'binomlog', inter_par2 = c(1, 20, 20))
funcMP_3D(as.numeric(m[1,]), distr = 'logbinom', inter_par2 = c(1, 20, 20))

# Проверка на стабильность модели
plot(x = 1:50, y = seq(0, 1, length.out = 50), col = "white")
for (j in 11:19){
  n <- 1:50
  res <- c()
  for (i in n){
    if (j <= 10){
      res <- c(res, hist_make(i, as.numeric(m[j,]), get_hist = FALSE, 'binomlog')[4])
    }
    else{
      res <- c(res, hist_make(i, as.numeric(m[j,]), get_hist = FALSE, 'logbinom')[4])
    }
  }
  lines(x = n, y = res, col = j)
}

# Корреляция с радиацией параметров и среднего
plot(seq(0, 45, 5), df$p[1:10], type = "l", col = "blue", ylim = c(0, max(df$p)), ylab = "Parametr p", xlab = "Dose, Gy", main = "Black - in vitro; Blue - in vivo")
lines(seq(0, 40, 5), df$p[11:19], type = "b", col = "black")

plot(seq(0, 45, 5), df$q[1:10], type = "l", col = "blue", ylim = c(0, max(df$q)), ylab = "Parametr q", xlab = "Dose, Gy", main = "Black - in vitro; Blue - in vivo")
lines(seq(0, 40, 5), df$q[11:19], type = "b", col = "black")

plot(seq(0, 45, 5), df$n[1:10], type = "l", col = "blue", ylim = c(0, max(df$n)), ylab = "Parametr n", xlab = "Dose, Gy", main = "Black - in vitro; Blue - in vivo")
lines(seq(0, 40, 5), df$n[11:19], type = "b", col = "black")

plot(seq(0, 45, 5), df$p[1:10] * df$n[1:10], type = "l", col = "blue", ylim = c(0, max(df$p * df$n)), ylab = "Average, n*p", xlab = "Dose, Gy", main = "Black - in vitro; Blue - in vivo")
lines(seq(0, 40, 5), df$p[11:19] * df$n[11:19], type = "b", col = "black")

mid <- c(0.38, 0.67, 0.83, 1.15, 1.70, 1.08, 1.83, 1.91, 2.23, 2.47, 0.39, 0.35, 0.59, 0.87, 0.57, 1.11, 1.13, 1.44, 1.69)

lines(seq(0, 45, 5), mid[1:10], type = "l", lty = 2, col = "green")
lines(seq(0, 40, 5), mid[11:19], type = "l", lty = 2, col = "black")

# Работа со словами
library(dplyr)
library(ggplot2)
library(foreach)
library(doSNOW)

WordSample <- read.csv("~/All/study_materials/Диплом/МоделированиеВR/WordCutter/output.csv")

num.word <- 1000

start_time <- Sys.time()

df.word.nbinom <- data.frame(name = row.names(WordSample)[1:num.word], n = rep(0, num.word), q = rep(0, num.word), p_value = rep(0, num.word))

for (i in 1:length(df.word.nbinom$name)){
  histWord <- hist(as.numeric(WordSample[i, ]), probability = TRUE, breaks = seq(-0.5, max(WordSample[i, ]) + 0.5, 1), xlim = c(-0.5, max(WordSample[i, ]) + 0.5), plot = FALSE)
  res <- hist_make(0, histWord$counts, get_hist = FALSE, distr = 'nbinom')
  df.word.nbinom[i,] <- c(df.word.nbinom[i, 1], res)
}

timediff <- difftime(Sys.time(),start_time)
cat("Расчёт занял: ", timediff, units(timediff))

df.word.nbinom <- data.frame()

start_time <- Sys.time()

cl <- makeCluster(5, type = "SOCK")
registerDoSNOW(cl)

df.word.nbinom <- foreach(i = 1:num.word, .combine = rbind, .inorder = TRUE) %dopar% {
  histWord <- hist(as.numeric(WordSample[i, ]), probability = TRUE, breaks = seq(-0.5, max(WordSample[i, ]) + 0.5, 1), xlim = c(-0.5, max(WordSample[i, ]) + 0.5), plot = FALSE)
  res <- hist_make(0, histWord$counts, get_hist = FALSE, distr = 'nbinom')
  c(row.names(WordSample)[i], res)
}

colnames(df.word.nbinom) <- c("name", "n", "q", "p_value")
df.word.nbinom <- as.data.frame(df.word.nbinom)

timediff <- difftime(Sys.time(),start_time)
cat("Расчёт занял: ", timediff, units(timediff))

df.word.nbinom <- df.word.nbinom |> mutate(n = as.numeric(n), q = as.numeric(q), p_value = as.numeric(p_value))

df.word.nbinom |> filter(n < 20) |> mutate(p_value_group = as.factor(ifelse(p_value < 0.1, 1, 2))) |>
  ggplot(aes(x = q, y = n, color = p_value_group)) +
  geom_point()

start_time <- Sys.time()

df.word.binlog <- data.frame(name = row.names(WordSample)[1:num.word], n = rep(0, num.word), q = rep(0, num.word), p = rep(0, num.word), p_value = rep(0, num.word))

for (i in 1:num.word){
  histWord <- hist(as.numeric(WordSample[i, ]), probability = TRUE, breaks = seq(-0.5, max(WordSample[i, ]) + 0.5, 1), xlim = c(-0.5, max(WordSample[i, ]) + 0.5), plot = FALSE)
  res <- hist_make(0, histWord$counts, get_hist = FALSE, distr = 'binomlog')
  df.word.binlog[i,] <- c(df.word.binlog[i, 1], res)
  print(i)
}

timediff <- difftime(Sys.time(),start_time)
cat("Расчёт занял: ", timediff, units(timediff))

df.word.binlog <- data.frame()

start_time <- Sys.time()

cl <- makeCluster(5, type = "SOCK")
registerDoSNOW(cl)

df.word.binlog <- foreach(i = 1:num.word, .combine = rbind, .inorder = TRUE) %dopar% {
  histWord <- hist(as.numeric(WordSample[i, ]), probability = TRUE, breaks = seq(-0.5, max(WordSample[i, ]) + 0.5, 1), xlim = c(-0.5, max(WordSample[i, ]) + 0.5), plot = FALSE)
  res <- hist_make(0, histWord$counts, get_hist = FALSE, distr = 'binomlog')
  c(row.names(WordSample)[i], res)
}

colnames(df.word.binlog) <- c("name", "n", "q", "p", "p_value")
df.word.binlog <- as.data.frame(df.word.binlog)

timediff <- difftime(Sys.time(),start_time)
cat("Расчёт занял: ", timediff, units(timediff))

name <- 'back'
word <- as.numeric(WordSample[name, ])
histWord <- hist(as.numeric(word), probability = TRUE, breaks = seq(-0.5, max(word) + 0.5, 1), xlim = c(-0.5, max(word) + 0.5), main = "", xlab = "Sample")
res <- hist_make(0, histWord$counts, get_hist = TRUE, distr = 'binomlog')
res
res <- hist_make(0, histWord$counts, get_hist = TRUE, distr = 'doublenbinom')
res
res <- hist_make(0, histWord$counts, get_hist = TRUE, distr = 'nbinom')
res

df.word.binlog <- df.word.binlog |> mutate(n = as.numeric(n), q = as.numeric(q), p = as.numeric(p), p_value = as.numeric(p_value))

df.word.binlog.better <- df.word.binlog[(df.word.nbinom$p_value > 0.2), ]
df.word.binlog.better <- df.word.binlog[(df.word.binlog$p_value > 0.2), ]
df.word.binlog.better <- df.word.binlog[(df.word.binlog$p_value > 0.2) & (df.word.nbinom$p_value > 0.2), ]

df.word.nbinom |> mutate(p_value_group = as.factor(ifelse((df.word.nbinom$p_value < df.word.binlog$p_value), 1, 2))) |> filter(n < 20) |>
  ggplot(aes(x = q, y = n, color = p_value_group)) +
  geom_point()

open3d()
lines3d(c(0, 1), c(0, 0), c(0, 0), color = "red")
lines3d(c(0, 0), c(0, 1), c(0, 0), color = "red")
points3d(df.word.binlog.better$q, df.word.binlog.better$p, log(df.word.binlog.better$n) / 10, color ="blue") 
points3d(df.word.binlog$q, df.word.binlog$p, log(df.word.binlog$n) / 10, color ="black") 

BadWord <- merge(df.word.nbinom |> select(-n, -q, -e), WordSample[1:num.word, ], by.x = "name", by.y = "word")
BadWord <- BadWord |> filter(p_value < 0.1)
write.csv(BadWord, file = "./BadWord.csv")

# Проверка оценки для свёртки отрицательных биномов
samDNB <- rdoublenbinom(1000, 0.8, 3, 0.5, 10)
histDNB <- hist(as.numeric(samDNB), probability = TRUE, breaks = seq(-0.5, max(samDNB) + 0.5, 1), xlim = c(-0.5, max(samDNB) + 0.5), plot = TRUE)
res <- hist_make(0, as.numeric(histDNB$counts), get_hist = TRUE, 'doublenbinom')
res <- hist_make(0, as.numeric(histDNB$counts), get_hist = TRUE, 'binomlog')
funcMP_3D(histDNB$counts, distr = 'binomlog')
res

