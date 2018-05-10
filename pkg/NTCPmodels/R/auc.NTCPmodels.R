auc.NTCPmodels<-function (y_pred, y_true)
{
  rank <- rank(y_pred)
  n_pos <- sum(y_true == 1)
  n_neg <- sum(y_true == 0)
  AUC <- (sum(rank[y_true == 1]) - n_pos * (n_pos + 1)/2)/(n_pos *
                                                             n_neg)
  return(AUC)
}