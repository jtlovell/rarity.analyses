spatial.genetic.graph<-function(pop.obj=gs,
                                subset.column="species",
                                tosubset="stricta",
                                stratum="pop",
                                spatial.data=ll,
                                node.size="n",
                                node.size.expander=7,
                                node.size.mean=0,
                                popgraph.alpha=0.2,
                                map.zoom=6
){
  if(is.null(subset.column) | is.null(tosubset)){
    subdata<-pop.obj
  }else{
    subdata<-pop.obj[pop.obj[,subset.column]==tosubset,]
  }
  data <- to_mv(subdata)
  pops <- subdata[,stratum]
  graph <- popgraph(x = data, groups = pops, alpha=popgraph.alpha)
  eb <- edge.betweenness.community(graph, modularity=T, directed=F)
  cat("\n",tosubset, ": modularity = ", max(eb[[5]]), " n.groups = ", length(unique(eb$membership)), sep="")
  member <- eb$membership
  color<-rainbow(length(unique(member)))[member]
  color[is.na(color)]<-"black"
  V(graph)$color <- color
  if(is.null(subset.column) | is.null(tosubset)){
    ll.sub<-ll
  }else{
    ll.sub<-ll[ll$sp==tosubset,]
  }
  ll.sub<-ll.sub[ll.sub$pop %in% V(graph)$name,]
  graph2 <- decorate_graph( graph, ll.sub, stratum="pop")
  location <- c(mean(V(graph2)$long), mean(V(graph2)$lat))
  map <- get_map(location,maptype="satellite", zoom=map.zoom)
  #get Ae for each population so to plot
  distHaversine <- function(long, lat){
    dlong = (long[2] - long[1])*pi/180
    dlat  = (lat[2] - lat[1])*pi/180
    R = 6371;
    a = sin(dlat/2)*sin(dlat/2) + cos(lat[1])*cos(lat[2])*sin(dlong/2)*sin(dlong/2)
    c = 2 * atan2( sqrt(a), sqrt(1-a) )
    d = R * c
    return(d) # in km
  }
  bb <- attr(map,"bb")
  sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                     lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                     lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                     lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))
  sbar$distance = distHaversine(long = c(sbar$lon.start,sbar$lon.end),
                                lat = c(sbar$lat.start,sbar$lat.end))
  ptspermm <- 2.83464567 
  wt<-E(graph2)$weight/max(E(graph2)$weight)
  p<-  ggmap( map ) + coord_map()+
    geom_edgeset( aes(x=long,y=lat), graph2, color="white", size=1, alpha=wt)+
    geom_nodeset( aes(x=long, y=lat, color=color, fill=color, size=size), graph2) +
    xlab("Longitude") + ylab("Latitude") +
    geom_segment(data = sbar,
                 aes(x = lon.start,
                     xend = lon.end,
                     y = lat.start,
                     yend = lat.end))+
    geom_text(data = sbar,
              aes(x = lon.start,
                  y = lat.start,
                  label = paste(format(distance, 
                                       digits = 4,
                                       nsmall = 2),
                                'km')),
              hjust = -.5,
              vjust = -.5,
              size = 8/ptspermm)
  print(p)
  return(p)
} 