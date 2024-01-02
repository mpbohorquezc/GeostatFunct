library(sf)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(rgdal)

map_variable_mexico <- function(datos, variable) {

   map <- read_sf("Data_frame/Mexico/AreaMetropolitana/mpiosutm.shp")

   plot1 <- ggplot() +
            geom_sf(data = map, size = 0.3) +
            geom_point(data = datos, aes_string(x = "X", y = "Y",
                                                color = variable)) +
            geom_text_repel(data = datos, aes_string(x = "X", y = "Y",
                                                label = "Estacion")) +
            scale_color_viridis_c(direction = -1) +
            labs(
                x = "",
                y = "",
                title = paste(variable, "map CDMX")
            ) +
            theme_void() +
            theme(
                plot.title = element_text(size = 14,
                                         face = "bold",
                                         colour = "black"),
                plot.margin = unit(c(1, 1, 1.5, 1.2), "cm"))

    return(plot1)

}

simple_scatter_plot <- function(daticos, variable1, variable2) {

    plot1 <- ggplot(as.data.frame(daticos),
                    aes_string(x = variable1, y = variable2,
                               color = as.factor(1))) +
                    geom_point() +
                    scale_colour_viridis_d() +
                    labs(
                        x = variable1,
                        y = variable2
                        ) +
                    theme_light() +
                    theme(legend.position = "none")

    return(plot1)

}

simple_hist_plot <- function(daticos, variable) {

    plot1 <- ggplot(as.data.frame(daticos),
                aes_string(x = variable, fill = as.factor(1))) +
                geom_density(alpha = 0.6) +
                # geom_histogram(aes(y = ..density..), alpha = 0.3) +
                scale_fill_viridis_d() +
                labs(
                    x = variable
                    ) +
                theme_light() +
                theme(legend.position = "none")

    return(plot1)

}

modelo_ideal_plot <- function(daticos, variabler) {

    plot <- subplot(simple_scatter_plot(daticos, "X", variabler),
                      simple_scatter_plot(daticos, "Y", variabler),
                      simple_scatter_plot(daticos, "altitud", variabler),
                      simple_hist_plot(daticos, variabler),
                      nrows = 2)
    return(plot)
}


# simple_hist_plot(datos, "PM10")
# simple_scatter_plot(datos, "Y", "PM10")
# simple_scatter_plot(datos, "altitud", "PM10")
# simple_scatter_plot(datos, "X", "PM10")

save_variable_map <- function(datos, variable) {

    plot <- map_variable_mexico(datos, variable)

    ggsave(filename = paste(variable, "map_CDMX.png", sep = ""),
           plot = plot,
           path = "Plots/Maps/")

}


prediction_plot <- function(g_object, variable) {

    map <- readOGR("Data_frame/Mexico/AreaMetropolitana/mpiosutm.shp")
    new <- sp::spsample(map, n = 10000, type = "regular")
    coordinates(new) ~ x1 + x2
    colnames(new@coords) <- c("X", "Y")

    predic <- predict(g_object, newdata = new)

    prediction <- data.frame(predic)

    pred <- paste(variable, ".pred", sep = "")

    plot1 <- plot_ly(
                 x = prediction[["X"]],
                 y = prediction[["Y"]],
                 z = prediction[[pred]],
                 type = "heatmap",
                 colorbar = list(title = "Prediction"),
                 reversescale = T

            ) %>%
            layout(
                   title = paste(variable, "prediction"),
                   xaxis = list(showticklabels = FALSE),
                   yaxis = list(showticklabels = FALSE),
                   scene = list(aspectration = list(x = 1, y = 1))
                   )
    return(plot1)

}

variance_plot <- function(g_object, variable) {

    map <- readOGR("Data_frame/Mexico/AreaMetropolitana/mpiosutm.shp")
    new <- sp::spsample(map, n = 10000, type = "regular")
    coordinates(new) ~ x1 + x2
    colnames(new@coords) <- c("X", "Y")

    predic <- predict(g_object, newdata = new)

    prediction <- data.frame(predic)

    var <- paste(variable, ".var", sep = "")

    plot2 <- plot_ly(
                 x = prediction[["X"]],
                 y = prediction[["Y"]],
                 z = prediction[[var]],
                 type = "heatmap",
                 colorbar = list(title = "variance"),
                 colorscale = "Hot",
                 reversescale = T
            )  %>%
            layout(
                   xaxis = list(showticklabels = FALSE),
                   yaxis = list(showticklabels = FALSE),
                   scene = list(aspectration = list(x = 1, y = 1))
                   )
    return(plot2)

}

cv_plot <- function(g_object, variable) {

    map <- readOGR("Data_frame/Mexico/AreaMetropolitana/mpiosutm.shp")
    new <- sp::spsample(map, n = 10000, type = "regular")
    coordinates(new) ~ x1 + x2
    colnames(new@coords) <- c("X", "Y")

    predic <- predict(g_object, newdata = new)

    prediction <- data.frame(predic)
    pred <- paste(variable, ".pred", sep = "")
    var <- paste(variable, ".var", sep = "")

    aux <- abs(sqrt(prediction[var]) / abs(prediction[pred]))
    aux[aux > 1] <- 1
    prediction["cv"] <- aux

    plot3 <- plot_ly(
                 x = prediction[["X"]],
                 y = prediction[["Y"]],
                 z = prediction[["cv"]],
                 type = "heatmap",
                 colorbar = list(title = "cv"),
                 colorscale = "Electric",
                 reversescale = T,
                 zmid = 0.5
            )  %>%
            layout(
                   xaxis = list(showticklabels = FALSE),
                   yaxis = list(showticklabels = FALSE),
                   scene = list(aspectration = list(x = 1, y = 1)))


    return(plot3)

}


mse <- function(g_object, variable, datos1) {

    datos <- datos1
    coordinates(datos) <- ~ X + Y
    predic <- predict(g_object, newdata = datos)

    prediction <- data.frame(predic)
    variable_a <- paste(variable, ".pred", sep = "")
    mse <- mean((prediction[[variable_a]] - datos[[variable]]) ** 2, na.rm = T)
    return(mse)

}