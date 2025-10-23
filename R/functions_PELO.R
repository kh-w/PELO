#' Calculate the winning probability of player 1
#'
#' @param p1_r Rating of player 1
#' @param p2_r Rating of player 2
#' @param p1_pw Pairwise advantage of of player 1 over player 2
#' @return The winning probability of player 1
#' @export
p_1win2___KK2 <- function(p1_r, p2_r, p1_pw, draw){return((10^(p1_r/400))/((10^(p1_r/400))+draw+(10^((p2_r-p1_pw)/400))))}

rating_pairwise_update___KK2 <- function(p1_r, p2_r, p1_pw, p2_pw, gameresult, K, K2, draw){

  p1_r_new <- p1_r
  p2_r_new <- p2_r
  p1_pw_new <- p1_pw
  p2_pw_new <- p2_pw
  p1win_prob <- p_1win2___KK2(p1_r, p2_r, p1_pw, draw)
  p2win_prob <- p_1win2___KK2(p2_r, p1_r, p2_pw, draw)
  draw_prob <- 1 - p1win_prob - p2win_prob

  game_result <- if(gameresult == "a"){1}else if(gameresult == "b"){0}else{0.5}

  p1_r_new <- p1_r_new + K*(game_result - p1win_prob)
  p2_r_new <- p2_r_new - K*(game_result - p1win_prob)
  p1_pw_new <- p1_pw_new + K2*(game_result - p1win_prob)
  p2_pw_new <- p2_pw_new - K2*(game_result - p1win_prob)

  return(c(p1_r_new, p2_r_new, p1_pw_new, p2_pw_new))

}

get_negloglik___KK2 <- function(KK2d, elo_type, data, elo_rat_KK2, output){

  df <- data
  elo_rating_KK2 <- elo_rat_KK2
  K <- KK2d[1]
  K2 <- KK2d[2]*if(elo_type=="elo"){0}else{if(elo_type=="pelo"){1}else{stop("specify elo or pelo")}}
  draw <- KK2d[3]*0

  for (row in 1:nrow(df)){

    player1 <- df$A[row]
    player2 <- df$B[row]

    player1_rating <- elo_rating_KK2$elo_rating_KK2[elo_rating_KK2$Player == player1]
    player2_rating <- elo_rating_KK2$elo_rating_KK2[elo_rating_KK2$Player == player2]

    player1_pairwise <- elo_rating_KK2[elo_rating_KK2$Player == player1, player2]
    player2_pairwise <- elo_rating_KK2[elo_rating_KK2$Player == player2, player1]

    player1_gplayed <- elo_rating_KK2$Games_played[elo_rating_KK2$Player == player1]
    player2_gplayed <- elo_rating_KK2$Games_played[elo_rating_KK2$Player == player2]

    if(df$match_result[row] == "a"){
      probability <- p_1win2___KK2(player1_rating, player2_rating, player1_pairwise, draw)
      probability <- max(0.000001,min(0.9999, probability))
      game_result <- df$match_result[row]
    }else if(df$match_result[row] == "b"){
      probability <- p_1win2___KK2(player2_rating, player1_rating, player2_pairwise, draw)
      probability <- max(0.000001,min(0.9999, probability))
      game_result <- df$match_result[row]
    }else{
      probability <- 1 - p_1win2___KK2(player2_rating, player1_rating, player2_pairwise, draw) - p_1win2___KK2(player1_rating, player2_rating, player1_pairwise, draw)
      probability <- max(0.000001,min(0.9999, probability))
      game_result <- df$match_result[row]
    }

    df$match_loglikelihood[row] <- log(probability)

    new_rating <- rating_pairwise_update___KK2(player1_rating, player2_rating, player1_pairwise, player2_pairwise, game_result, K, K2, draw)

    elo_rating_KK2$elo_rating_KK2[elo_rating_KK2$Player == player1] <- new_rating[1]
    elo_rating_KK2$elo_rating_KK2[elo_rating_KK2$Player == player2] <- new_rating[2]
    elo_rating_KK2[elo_rating_KK2$Player == player1, player2] <- new_rating[3]
    elo_rating_KK2[elo_rating_KK2$Player == player2, player1] <- new_rating[4]

    elo_rating_KK2$Games_played[elo_rating_KK2$Player == player1] <- player1_gplayed + 1
    elo_rating_KK2$Games_played[elo_rating_KK2$Player == player2] <- player2_gplayed + 1

  }

  if(output == "negloglik"){
    print(c(-sum(df$match_loglikelihood), KK2d))
    return(-sum(df$match_loglikelihood))
  } else if(output == "all"){
    df$match_likelihood <- df$match_loglikelihood %>% exp()
    return(list("df" = elo_rating_KK2, "negloglik" = -sum(df$match_loglikelihood), "prob" = df))
  } else {
    stop("specify output as 'negloglik' or 'all'.")
  }

}

prediction <- function(p1_rating, p2_rating, p1_pairwise_p2, draw, prob = FALSE){
  p <- p_1win2___KK2(p1_r = p1_rating, p2_r = p2_rating, p1_pw = p1_pairwise_p2, draw = draw)
  q <- p_1win2___KK2(p1_r = p2_rating, p2_r = p1_rating, p1_pw = -p1_pairwise_p2, draw = draw)
  d <- 1 - p - q
  if(prob == TRUE){
    return(max(p,q,d))
  }else{
    if(max(p,q,d) == p){
      return("a")
    }else if(max(p,q,d) == q){
      return("b")
    }else{
      return("d")
    }
  }
}

get_optimized_para <- function(initial_guess, elo_type, data, elo_rat_KK2){
  result <- optim(initial_guess,
                  fn = get_negloglik___KK2,
                  method = "L-BFGS-B",
                  lower = c(0.01,-10^7,0.01),
                  upper = c(10^7,10^7,10^7),
                  control = list(trace = 1),
                  elo_type = elo_type,
                  data = data,
                  elo_rat_KK2 = elo_rat_KK2,
                  output = "negloglik")$par
  return(result)
}

likelihood_ratio_test <- function(initial_guess, data, elo_rat_KK2){
  tic <- Sys.time()

  optimized_para_KK2_elo <- get_optimized_para(initial_guess = initial_guess,
                                               elo_type = "elo",
                                               data = data,
                                               elo_rat_KK2 = elo_rat_KK2)
  optimized_para_KK2_pelo <- get_optimized_para(initial_guess = initial_guess,
                                                elo_type = "pelo",
                                                data = data,
                                                elo_rat_KK2 = elo_rat_KK2)

  toc <- Sys.time()

  runtime_KK2 <- toc - tic
  print(runtime_KK2)

  negloglik_optimized_KK2_pelo <- get_negloglik___KK2(KK2d = optimized_para_KK2_pelo,
                                                      elo_type = "pelo",
                                                      data = data,
                                                      elo_rat_KK2 = elo_rat_KK2,
                                                      output = "all")
  negloglik_optimized_KK2_elo <- get_negloglik___KK2(KK2d = optimized_para_KK2_elo,
                                                      elo_type = "elo",
                                                      data = data,
                                                      elo_rat_KK2 = elo_rat_KK2,
                                                      output = "all")

  lrt_statistic <- -2*(negloglik_optimized_KK2_pelo$negloglik - negloglik_optimized_KK2_elo$negloglik)
  return(c("K_mle" = optimized_para_KK2_elo, "KK2_mle" = optimized_para_KK2_pelo, "nll_elo" = negloglik_optimized_KK2_elo$negloglik, "nll_pelo" = negloglik_optimized_KK2_pelo$negloglik, "lrt_statistic" = lrt_statistic))
}
