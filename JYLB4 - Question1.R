library(here)
library(sf)
library(tmap)
library(dplyr)
library(units)
library(sp)
library(spdep)
library(reshape2)
library(spatialreg)


#load data
TL_greenspace <- st_read(here::here("Desktop/GEOG0114/greenspace", "TL_GreenspaceSite.shp"))
TQ_greenspace <- st_read(here::here("Desktop/GEOG0114/greenspace", "TQ_GreenspaceSite.shp"))
ward <- st_read(here::here("Desktop/GEOG0114/statistical-gis-boundaries-london/ESRI", "London_Ward_CityMerged.shp"))
LSOA <- st_read(here::here("Desktop/GEOG0114/statistical-gis-boundaries-london/ESRI","LSOA_2011_London_gen_MHW.shp"))
LSOA <- LSOA[,-c(3:8)]
ward <- st_buffer(ward, 0)
london_outline <- ward %>% summarise(area = sum(HECTARES))


#select greenspace within London
london_greenspace = all_greenspace[london_outline,] %>% st_intersection(london_outline)

#caculate area of each greenspace
london_greenspace$area_m <- st_area(london_greenspace)

#create four sizes of buffers
ha2_gs_buffer_300m <- st_buffer(london_greenspace %>% filter(area_m > set_units(20000.0, m^2)), dist = 300)
ha20_gs_buffer_2km <- st_buffer(london_greenspace %>% filter(area_m > set_units(200000.0, m^2)), dist = 2000)
ha100_gs_buffer_5km <- st_buffer(london_greenspace %>% filter(area_m > set_units(1000000.0, m^2)), dist = 5000)
ha500_gs_buffer_10km <- st_buffer(london_greenspace %>% filter(area_m > set_units(5000000.0, m^2)), dist = 10000)

ha2_gs_buffer_300m_single <- ha2_gs_buffer_300m %>% summarise(area = sum(st_area(ha2_gs_buffer_300m)))
ha20_gs_buffer_2km_single <- ha20_gs_buffer_2km %>% summarise(area = sum(st_area(ha20_gs_buffer_2km)))
ha100_gs_buffer_5km_single <- ha100_gs_buffer_5km %>% summarise(area = sum(st_area(ha100_gs_buffer_5km)))
ha500_gs_buffer_10km_single <- ha500_gs_buffer_10km %>% summarise(area = sum(st_area(ha500_gs_buffer_10km)))

#generate four series of boolean variables to measure greenspace accessibility
LSOA$ha2_access <- st_intersects(lsoa_centorid, ha2_gs_buffer_300m_single, sparse=FALSE)
LSOA$ha20_access <- st_intersects(lsoa_centorid, ha20_gs_buffer_2km_single, sparse=FALSE)
LSOA$ha100_access <- st_intersects(lsoa_centorid, ha100_gs_buffer_5km_single, sparse=FALSE)
LSOA$ha500_access <- st_intersects(lsoa_centorid, ha500_gs_buffer_10km_single, sparse=FALSE)

#plot greenspace accessibility in each LSOA
tm_shape(LSOA) + tm_fill(col = 'ha2_access',palette = 'Greens')+tm_layout(frame =  F)
tm_shape(LSOA) + tm_fill(col = 'ha20_access',palette = 'Greens')+tm_layout(frame =  F)
tm_shape(LSOA) + tm_fill(col = 'ha100_access',palette = 'Greens')+tm_layout(frame =  F)
tm_shape(LSOA) + tm_fill(col = 'ha500_access',palette = 'Greens')+tm_layout(frame =  F)

#genarate MSOA-scale data 
LSOA$ha2_access <-  as.integer(as.logical(LSOA$ha2_access))
LSOA$ha20_access <-  as.integer(as.logical(LSOA$ha20_access))
LSOA$ha100_access <-  as.integer(as.logical(LSOA$ha100_access))

LSOA <- st_buffer(LSOA, 0)

LSOA$tmp <- rep(1,nrow(LSOA))
MSOA_access <- LSOA %>% group_by(MSOA11CD) %>% summarise(across(c(ha2_access,ha20_access,ha100_access,tmp), sum))
MSOA_access$avg_acess <- (MSOA_access$ha2_access + MSOA_access$ha20_access + MSOA_access$ha100_access) / MSOA_access$tmp

#plot AAS
tm_shape(MSOA_access)+tm_fill(col = 'avg_acess')


#GLOBAL moran'I
#use k-nearest neighbours of 4

coordsM <- MSOA_access %>% st_centroid()%>% st_geometry()
knn_MSOA <-coordsM %>% knearneigh(., k=4)
LMSOA_knn <- knn_MSOA %>% knn2nb()
LMSOA.knn_4_weight <- LMSOA_knn %>% nb2listw(., style="C")

moran.test(MSOA_access$avg_acess,LMSOA.knn_4_weight)


#Demogrphic traits
msoa_profile <- read.csv('Desktop/GEOG0114/msoa-data.csv',check.names = F)
head(msoa_profile)

MSOA_access <- left_join(MSOA_access,msoa_profile,by = c('MSOA11CD' = 'Middle Super Output Area'))


#correlation test for Age 15-
cor.test(MSOA_access$`Age0-15`,MSOA_access$avg_acess)

#correlation test for Age16-44
cor.test(MSOA_access$`Age16-29`,MSOA_access$avg_acess)
cor.test(MSOA_access$`Age30-44`,MSOA_access$avg_acess)

#correlation test for Age45+
cor.test(MSOA_access$`Age45-64`,MSOA_access$avg_acess)
cor.test(MSOA_access$`Age65+`,MSOA_access$avg_acess)

#scatterplot
MSOA_plot <- cbind(MSOA_access$`Age0-15`,MSOA_access$`Age16-29`+MSOA_access$`Age30-44`,MSOA_access$`Age45-64`+MSOA_access$`Age65+`,MSOA_access$avg_acess)
MSOA_plot <- as.data.frame(MSOA_plot)
colnames(MSOA_plot) <- c('Age0-15','Age16-44','Age45+','AAS')

MSOA_plot <- melt(MSOA_plot, id.vars="AAS")
ggplot(MSOA_plot, aes(value,AAS, col=variable)) + geom_point() + scale_color_brewer(palette="Accent")+ labs( x = 'population proportion')



#spatially-lagged regression model

obesity_model <- lagsarlm(MSOA_access$`Obesity population aged 16+` ~ MSOA_access$avg_acess + MSOA_access$`Day-to-day activities not limited` + MSOA_access$white_proportion, data = MSOA_access,
                          LMSOA.knn_4_weight,
                          method = "eigen")

summary(obesity_model)

MSOA_access$residuals <- residuals(obesity_model)
moran.test(MSOA_access$residuals,LMSOA.knn_4_weight)
