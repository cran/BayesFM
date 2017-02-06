#' @export

print.befa <- function(x, ...)
{
  if (class(x) != 'befa')
    stop('object passed to print.befa should be of class befa')

  cat('BEFA output - Bayesian Exploratory Factor Analysis\n\n')
  cat('call:\n  ')
  print(attr(x, 'call'))
  cat('\n')
  cat('posterior column switch:', attr(x, "post.column.switch"), '\n')
  cat('posterior sign switch:  ', attr(x, "post.sign.switch"), '\n\n')
}
