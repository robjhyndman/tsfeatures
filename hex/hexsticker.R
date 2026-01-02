library(cropcircles)
library(magick)
library(tibble)
library(ggplot2)
library(ggpath)
library(showtext)

# choose a font from Google Fonts
font_add_google("Fira Sans", "firasans")
showtext_auto()

img_padded <- image_border(
  image_read(here::here("hex/pcs.png")),
  color    = "transparent",
  geometry = "140x200" 
) |>
image_colorize(opacity = 35, color = "white")


# Write padded image to a temporary file (or a real file if you prefer)
padded_path <- here::here("hex/pcs_padded.png")
image_write(img_padded, padded_path)

img_cropped <- hex_crop(
  images = padded_path,
  bg_fill = "#f5f4f3ff",
  border_colour = "#c14b14",
  border_size = 50
)

ggplot() +
  geom_from_path(aes(0.5, 0.5, path = img_cropped)) +
  annotate(
    "text",
    x = 0.20,
    y = 0.23,
    label = "tsfeatures",
    family = "firasans",
    size = 36,
    colour = "#c14b14",
    angle = 0,
    hjust = 0,
    fontface = "bold"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  coord_fixed()

ggsave("./man/figures/tsfeatures-hex.png", height = 3, width = 3)

# Trim transparent edges
img <- image_read("./man/figures/tsfeatures-hex.png")
img_trim <- image_trim(img)

image_write(img_trim, "./man/figures/tsfeatures-hex.png")
