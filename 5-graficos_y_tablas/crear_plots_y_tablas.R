pacman::p_load(vcfR, dplyr, stringr, ggplot2)

tuarchivo <- "muestra1.strelka.variants_VEP.ann.vcf.gz"

vcf <- read.vcfR( tuarchivo )

# Extract the INFO column
info <- vcf@fix[, "INFO"]

# Extract CSQ field from INFO
csq_data <- str_extract(info, "CSQ=[^;]+") %>%
  str_remove("^CSQ=") %>%
  str_split_fixed(",", n = 1)  # Get only the first annotation per site

# Split CSQ string by '|'
csq_split <- str_split_fixed(csq_data, "\\|", n = 13)
csq_df <- as.data.frame(csq_split)
colnames(csq_df)[2] <- "Consequence"  # 2nd field is Consequence

# Filter out missing values
csq_df <- csq_df %>% filter(Consequence != "")

# Classify variant effects
csq_df <- csq_df %>%
  mutate(
    type = case_when(
      str_detect(Consequence, "missense|stop_gained|frameshift|start_lost|stop_lost|inframe|protein_altering") ~ "Non-synonymous",
      str_detect(Consequence, "synonymous_variant") ~ "Synonymous",
      str_detect(Consequence, "intron|intergenic|upstream|downstream|non_coding|UTR") ~ "Non-coding",
      TRUE ~ "Other"
    ),
    broad_category = case_when(
      type %in% c("Synonymous", "Non-synonymous") ~ "Coding",
      type == "Non-coding" ~ "Non-coding",
      TRUE ~ "Other"
    )
  )

# Count for plotting
plot_df <- csq_df %>%
  count(broad_category, type) %>%
  group_by(broad_category) %>%
  mutate(
    fraction = n / sum(n),
    label = paste0(type, " (", n, ")")
  )

# Donut pie plot
ggplot(plot_df, aes(x = 2, y = fraction, fill = type)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ broad_category) +
  xlim(0.5, 2.5) +  # creates donut shape
  theme_void() +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Variant Effects: Coding vs Non-coding") +
  theme(legend.position = "none")

#####
# Load required package
library(circlize)
library(dplyr)
library(stringr)

# Ensure your data looks like this:
# Cargar paquetes
pacman::p_load( "vcfR",        # Para leer el vcf
                "dplyr",       # Para filtrar data
                "stringr",     # Para detectar palabras
                "circlize" )   # Para dibujar el circos

# leer las variantes
vcf_file <- vcf         # leemos directamente el vcf desde internet

# extraemos la informacion de cromosoma y posicion de todas las variantes
variant.df <- data.frame(
  chr = vcf_file@fix[ , "CHROM" ],              # el objecto vcf_file tiene cada columna del VCF, usamos la CHROM
  start = vcf_file@fix[ , "POS" ] %>%             # y la columna POS
    as.numeric( )
) %>% 
  mutate( variant_val = 1 )


# Preparar los parametros del Circos
circos.par( start.degree = 70,                         # controlamos el angulo en el que comienza el circos
            gap.degree = c( rep( 1, 21 ), 40 ) )       # preparamos 21 espacios de 1 grado, y un ultimo espacio de 40 grados

# inicializamos el circos
cromosomas_circos <- paste0( "chr", 1:22 )             # prepara los que quieres dibujar

circos.initializeWithIdeogram( species = "hg38",                        # especie humana, version GRCh38
                               track.height = 0.05,                     # grosor del 5% en el track
                               chromosome.index = cromosomas_circos )   # aqui indicas que cromosomas quieres usar

# empezar el track para las variantes
# encontrar el valor maximo de densidad de variantes
max_valor <- max( variant.df$variant_val ) 

circos.track(
  
  ylim = c( 0, max_valor ) %>% rev( ),     # el rango del eje se invierte para lograr un efecto de dibujo inverso
  bg.border = NA,                          # borramos el borde de las celdas
  track.height = 0.3                       # el track ocupara un 30% del radio del circos
  
)

# iniciamos un loop for para dibujar cada cromosoma unico
for ( chrom in unique( variant.df$chr ) ) {
  
  # creamos un dataframe intermedio para cada cromosoma en turno
  chrom_data <- variant.df %>%
    filter( chr == chrom )                 # filtramos solo filas donde el cromosoma sea el crom en turno
  
  # dibujamos las variantes para el cromosoma en turno  
  circos.points(
    x = chrom_data$start,
    y = chrom_data$variant_val,
    # col = "skyblue",
    pch = 16,
    sector.index = chrom                   # indicamos que se dibuja el cromosoma en turno
  )
  
}

# Cerramos el circos
circos.clear()

#
# pass the VCF to excel ######
# Install if necessary
# install.packages("readr")
# install.packages("writexl")
# install.packages("stringr")
# install.packages("dplyr")

# Load packages
library(readr)
library(writexl)
library(stringr)
library(dplyr)

# Path to your file
vcf_file <- "muestra1.strelka.variants_VEP.ann.vcf.gz"  # or "yourfile.vcf.gz" if compressed

# Read lines, skipping metadata lines that start with ##
vcf_lines <- read_lines(vcf_file)
vcf_data <- vcf_lines[!str_starts(vcf_lines, "##")]

# Split the header and data
header_line <- vcf_data[1]
data_lines <- vcf_data[-1]

# Parse header
colnames <- str_split(header_line, "\t")[[1]]

# Split each line into fields
data_split <- str_split(data_lines, "\t")

# Convert to data frame
vcf_df <- data_split %>%
  do.call(rbind, .) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Set column names
colnames(vcf_df) <- colnames

write_xlsx(vcf_df, "parsed_vcf.xlsx")
