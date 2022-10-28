## Figure
# If we want to use this to make a figure, we can pass a sequence
# of values to x1. Remember, x2 must always be 1 row. Let's 
# visualize the RSS for forage vs. mean forage.
x1 <- data.frame(forage = seq(0, 800, length.out = 100), 
                 predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

x2 <- data.frame(forage = mean(values(hab$forage)), 
                 predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

logRSS <- log_rss(m1, x1, x2, ci = "se", ci_level = 0.95)

# We have a plot method for 'log_rss' objects to make a very basic figure.
plot(logRSS)

# But if we want more control, we can use ggplot with the 'df' data.frame.
ggplot(logRSS$df, aes(x = forage_x1, y = exp(log_rss), 
                      ymin = exp(lwr), ymax = exp(upr))) +
  geom_ribbon(color = "black", fill = "gray80", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  xlab("Forage at x1") +
  ylab("RSS vs. Mean Forage") +
  theme_bw()