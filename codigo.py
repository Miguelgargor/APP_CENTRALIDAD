### A) LIBRERÍAS NECESARIAS:
import pandas as pd                # Tratamiento de datos.
import matplotlib.pyplot as plt    # Ploteo de Gráficos.
import osmnx as ox                 # Mapas.
import networkx as nx              # Grafos.
from adjustText import adjust_text # Para ajustar las etiquetas de las aristas (Texto).
import os                          # Tratamiento de archivos y directorios.
import threading                   # Para poner el LÍMITE de TIEMPO para cargar el Grafo.
import geopy                       # Para obtener las COORDENADAS del LUGAR.
import geopandas as gpd            # Mapa interactivo.
import streamlit as st             # Creación de la APP WEB.
##=========================================================================================================================================##

## TÍTULO de la PESTAÑA del NAVEGADOR y añadir ICONO: (TÍTULO, URL ICONO (Subido a GitHub))
URL_ICONO= 'https://raw.githubusercontent.com/Miguelgargor/APP_CENTRALIDAD/main/Imagen1.png?token=GHSAT0AAAAAACFRVEI7FPDBRK67HLLDMTIKZGCREIA'
st.set_page_config(page_title="ANÁLISIS CENTRALIDAD", page_icon=URL_ICONO, layout="centered", initial_sidebar_state="auto")
##===========================================================================================================================================##
##===========================================================================================================================================##

#### FUNCIÓN PARA DIBUJAR EL GRAFO CON EL NODO DE MAYOR CENTRALIDAD ####
def Dibujar_Nodo_Central(Lugar, Tiempos_Viaje, Velocidad, Distancia=3500, Tiempo_Espera=180, Tipo_Calles='drive',
                         Tamaño_nodo=40, Tamaño_figura=(16,16), Guardar_img=True, Guardar_Grafo=True):
### B) CREACIÓN Y GUARDADO DEL GRAFO:
## 1º) Creación del Grafo:
    # 0) Si se pasa del tiempo saltará una Excepción:
    class TimeoutException(Exception):
        pass
    ##.................................................................................................................................##

    # 1) Carga del GRAFO desde un LUGAR:
    def lugar(Lugar, Tipo_Calles):
        def lugar_thread(result_list, Lugar, Tipo_Calles):
            try:
                g1 = ox.graph_from_place(Lugar,                    # Lugar del Grafo.
                                        truncate_by_edge=True,    # Sólo aristas dentro del Grafo.
                                        network_type=Tipo_Calles) # Tipo de Calles.
                result_list.append(g1) # Si se consigue cargar el Grafo dentro del TIEMPO LÍMITE, añádelo a la lista.

            except TimeoutException: # Si NO se consigue cargar el Grafo dentro del TIEMPO LÍMITE, devuelve la EXCEPCIÓN.
                print("Revasa el límite de tiempo")
        result_list= [] # Crear lista para añadir el Grafo.
        
        # Ejecución de la carga del GRAFO desde un LUGAR (con el límite de tiempo): [Pasándole todos los argumentos]:
        t= threading.Thread(target=lugar_thread, args=(result_list, Lugar, Tipo_Calles)); t.start()
        t.join(timeout=Tiempo_Espera)  # Espera a que el hilo termine ó exceda el límite de tiempo.
        
        # Si al cargar el GRAFO desde un LUGAR se pasa del límite de tiempo, sacar la Excepción:
        if t.is_alive():
        #  print("... Tiempo de espera excedido, cargando el Grafo desde el Punto ...")
            raise TimeoutException()
        # Si NO se pasa del LÍMITE de TIEMPO...
        if result_list:
            return result_list[0]  # ... Devuelve el Grafo.
        else:
            return None
    ##.................................................................................................................................##

    # 2a) Se cargará el GRAFO desde un LUGAR...
    try:
        g1= lugar(Lugar, Tipo_Calles)
    except TimeoutException:           # ... Y en caso de pasarse del LÍMITE de TIEMPO...
        g1= None                       # g1 = NoneType.
        print('... Tiempo de espera excedido, cargando el Grafo desde el Punto ...')

    # 2b) Se cargará el GRAFO desde un PUNTO:
    # 2b.1) Obtención de las COORDENADAS del LUGAR indicado:
    location= geopy.Nominatim(user_agent="my_app").geocode(Lugar)
    Coordenadas= (location.latitude, location.longitude)

    # 2b.2) GRAFO desde el PUNTO de las COORDENADAS:
    g2 = ox.graph_from_point(center_point=Coordenadas, # Coordenadas desde las que se creará el Grafo.
                            dist=Distancia,           # Distancia alrededor del punto de las Coordenadas anteriores.
                            dist_type="bbox",         # Creará una Caja con esa distancia.
                            network_type=Tipo_Calles) # Tipo de Calles.
    ##.................................................................................................................................##

    # 3) Selección del GRAFO que se usará en el resto del ESTUDIO:
    if g1 is not None and g2 is not None:
        if len(g1.nodes)+len(g1.edges) <= len(g2.nodes)+len(g2.edges):
            G= g1
            print('*** El Grafo que se usará en el estudio es el del LUGAR indicado, al haberse cargado dentro del límite de tiempo y ser el más compacto.')
        else:
            G= g2
            print('*** El Grafo que se usará en el estudio es el del PUNTO de COORDENADAS indicadas, al haberse cargado dentro del límite de tiempo y ser el más compacto.')

    elif g1 is None and g2 is not None:
        G= g2
        print('*** El Grafo que se usará en el estudio es el del PUNTO de COORDENADAS indicadas, al sólo cargarse éste dentro del límite de tiempo.')
    ##.................................................................................................................................##

    G_proj= ox.project_graph(G)                                                                # Proyección Geográfica del Grafo.
    G= ox.consolidate_intersections(G_proj, rebuild_graph=True, tolerance=15, dead_ends=False) # Simplificar el Grafo.
##------------------------------------------------------------------------------------------------------------------##
## 2º) Comprobar si la carpeta no existe antes de crearla:
    Nombre_Lugar= Lugar.split(',')[0].strip().replace(' ', '_') + '_' + Tipo_Calles # Nombre del Lugar extraído del parámetro 'Lugar'.
    
    if not os.path.exists('./'): # Si no existe ya la carpeta...GENERAL.
        os.makedirs('./ANALISIS CENTRALIDAD') # Creala con el Nombre del Lugar a analizar.

    if not os.path.exists('./'): # Si no existe ya la carpeta...ESPECÍFICA.
        os.makedirs('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}'.format(Nombre_Lugar)) # Creala con el Nombre del Lugar a analizar.
##------------------------------------------------------------------------------------------------------------------##
## 3º) Guardar Grafo:
    if Guardar_Grafo==True:
        ox.save_graphml(G, './ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/Grafo_{}.graphml'.format(Nombre_Lugar, Nombre_Lugar))
##=========================================================================================================================================##
##=========================================================================================================================================##

### C) CENTRALIDADES Y NODO CON MAYOR CENTRALIDAD MEDIA:
## 4º) Centralidad de INTERMEDIACIÓN:
    bc= nx.betweenness_centrality(G, weight='length')
    nx.set_node_attributes(G, bc, 'Centralid_Intermedia')
##------------------------------------------------------------------------------------------------------------------##
## 5º) Centralidad de Grado de ENTRADA:
    cen_Entrada= nx.in_degree_centrality(G)
    nx.set_node_attributes(G, cen_Entrada, 'Centralid_Entrada')
##------------------------------------------------------------------------------------------------------------------##
## 6º) Centralidad de Grado de SALIDA:
    cen_Salida= nx.out_degree_centrality(G)
    nx.set_node_attributes(G, cen_Salida, 'Centralid_Salida')
##------------------------------------------------------------------------------------------------------------------##
## 7º) Centralidad de Cercanía:
    cen_cercanía= nx.closeness_centrality(G, distance='length')
    nx.set_node_attributes(G, cen_cercanía, 'Centralid_Cercanía')
##------------------------------------------------------------------------------------------------------------------##
## 8º) RESUMEN CENTRALIDAD y Nodo_Clave (MÁS CENTRALIDAD):
    Intermedia= pd.DataFrame(bc.items(), columns=['Nodo','Cen_Intermedia'])       # Centralidad de Intermediación.
    Entrada= pd.DataFrame(cen_Entrada.items(), columns=['Nodo','Cen_Entrada'])    # Centralidad de Grado de Entrada.
    Salida= pd.DataFrame(cen_Salida.items(), columns=['Nodo','Cen_Salida'])       # Centralidad de Grado de Salida.
    Cercanía= pd.DataFrame(cen_cercanía.items(), columns=['Nodo','Cen_Cercanía']) # Centralidad de Cercanía.
    # a) Resumen:
    RESUMEN_CENTRALIDAD= Intermedia.merge(Entrada, on='Nodo').merge(Salida, on='Nodo').merge(Cercanía, on='Nodo')
    # b) Centralidad Media:
    RESUMEN_CENTRALIDAD['Promedio']= RESUMEN_CENTRALIDAD[['Cen_Intermedia', 'Cen_Entrada', 'Cen_Salida', 'Cen_Cercanía']].mean(axis=1, numeric_only=True)
    # c) NODO DE MAYOR CENTRALIDAD MEDIA:
    maximo= RESUMEN_CENTRALIDAD.loc[RESUMEN_CENTRALIDAD['Promedio'].idxmax()]
    Nodo_Clave= maximo['Nodo']

    ## ADICIÓN COMO ATRIBUTO DE LOS NODOS LA CENTRALIDAD MEDIA:
    nx.set_node_attributes(G, RESUMEN_CENTRALIDAD['Promedio'], 'Media_Centralidad')
##=========================================================================================================================================##
##=========================================================================================================================================##

### D) GRAFO CON EL NODO_CLAVE RESALTADO:
## 9º) Aristas incidentes en el Nodo Clave:
    aristas_incidentes = []
    ## A) Para el nodo 'u':
    for u, aristas_vecinas in G.adj.items():
        if u==Nodo_Clave:
            for v, key_dict in aristas_vecinas.items():
                for clave, Datos in key_dict.items():
                    if 'name' in Datos:  # Saca el atributo 'name' (Nombre de la calle, si lo tiene).
                        aristas_incidentes.append((u, v, clave, Datos['name']))
    ## B) Para el nodo 'v':   
        else:
            for v, key_dict in aristas_vecinas.items():
                if v==Nodo_Clave:
                    for clave, Datos in key_dict.items():
                        if 'name' in Datos:  # Saca el atributo 'name' (Nombre de la calle, si lo tiene).
                            aristas_incidentes.append((u, v, clave, Datos['name']))
##------------------------------------------------------------------------------------------------------------------##
## 10º) DIBUJAR el grafo con el Nodo Clave en color ROJO:
    fig, ax = ox.plot_graph(G,                                              # Grafo.
                            figsize=Tamaño_figura,                          # Tamaño figura.
                            node_color=['r' if node==Nodo_Clave else 'skyblue' for node in G.nodes], # Color nodos.
                            node_size=[Tamaño_nodo if node==Nodo_Clave else 1 for node in G.nodes],  # Tamaño nodos.
                            show=False, close=False, bgcolor='whitesmoke')  # Color fondo.
    plt.title('Nodo con mayor Centralidad Media en {}'.format(Nombre_Lugar), fontweight='bold') # Título Subgrafo con Zoom.
##------------------------------------------------------------------------------------------------------------------##
## 11º) Etiquetas de las aristas incidentes en el Nodo clave:
    etiquetas_aristas= []
    for _, _, _, nombre in aristas_incidentes:
        if nombre not in etiquetas_aristas:      # Conjunto de etiquetas únicas.
            etiquetas_aristas.append(nombre)

    # a) Agregar a la Leyenda las etiquetas:
    for etiqueta in etiquetas_aristas:
        ax.plot([], [], label=etiqueta, color='red')
    # b) Leyenda:
    ax.legend(title='Calles del Nodo Central', loc='best', fontsize=6)
    st.pyplot(fig)  # Mostrar el ploteo en la APP WEB.

## GUARDAR IMAGEN:
    if Guardar_img==True:
            plt.savefig('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/Grafo_{}_Nodo_Mayor_Centralidad.png'.format(Nombre_Lugar, Nombre_Lugar))
##------------------------------------------------------------------------------------------------------------------##
## 12º) Crear MAPA INTERACTIVO con el NODO CENTRAL RESALTADO EN ROJO:
    # a) Convertir el Grafo a un GeoDataFrame:
    geomapa_nodos = ox.graph_to_gdfs(G, edges=False)

    # b) Resaltar el Nodo_Clave en ROJO (y el resto en azul):
    geomapa_nodos['color']= ['red' if node_id == Nodo_Clave else 'blue' for node_id in geomapa_nodos.index]

    # c) Resaltar el Nodo_Clave con MAYOR TAMAÑO (y el resto más pequeño):
    geomapa_nodos['geometry'] = geomapa_nodos['geometry'].buffer(15)
    geomapa_nodos.loc[geomapa_nodos.index == Nodo_Clave, 'geometry'] = geomapa_nodos.loc[geomapa_nodos.index == Nodo_Clave, 'geometry'].buffer(30)

    # d) Mostrar como información el Nº de Nodo:
    geomapa_nodos['Nº Node']= geomapa_nodos.index

    # e) MAPA INTERACTIVO:    (Color_Nodos, Leyenda,     Estilo_Borde_Nodos,                        Estilo del Mapa,          Información al Pulsar, Resto de Info=NO):
    mapa_nodos = geomapa_nodos.explore(color='color', legend=True, style_kwds={'color':'black','weight':0.3}, tiles='CartoDB positron', popup=['Nº Node'], tooltip=None)

    # f) Verificar si la carpeta 'MAPAS_INTERACTIVOS' ya existe dentro de 'ANALISIS_CENTRALIDAD':
    if not os.path.exists('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/MAPAS_INTERACTIVOS'.format(Nombre_Lugar)):
        os.makedirs('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/MAPAS_INTERACTIVOS'.format(Nombre_Lugar))  # Si no existe-> crearla.

    # g) Guardar el Mapa Interactivo como archivo HTML (Google):
    mapa_nodos.save('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/MAPAS_INTERACTIVOS/MAPA_INTERACTIVO_Nodo_Mayor_Centralidad.html'.format(Nombre_Lugar))
##=========================================================================================================================================##
##=========================================================================================================================================##

#### E) SUBGRAFO CON ZOOM AL NODO_CLAVE:
## 13º) CREAR: Subgrafo con Zoom sobre el Nodo_Clave:
    vecinos= list(G.neighbors(Nodo_Clave))          # Vecinos del Nodo_Clave.
    Subgrafo_Zoom= G.subgraph(vecinos+[Nodo_Clave]) # Subgrafo (Vecinos + Nodo_Clave).
##------------------------------------------------------------------------------------------------------------------##
## 14º) DIBUJAR: Subgrafo con el Nodo_Clave en color ROJO:
    fig, ax= ox.plot_graph(Subgrafo_Zoom,                                 # Subgrafo.
                           figsize=(Tamaño_figura[0]-6, Tamaño_figura[1]-6), # Tamaño figura.
                           node_color=['r' if node == Nodo_Clave else 'skyblue' for node in Subgrafo_Zoom.nodes], # Color nodos.
                           node_size=[Tamaño_nodo if node == Nodo_Clave else Tamaño_nodo-10 for node in Subgrafo_Zoom.nodes], # Tamaño nodos.
                           show=False, close=False, bgcolor='whitesmoke') # Color fondo.
    plt.title('Zoom sobre el Nodo con mayor Centralidad Media en {}'.format(Nombre_Lugar), fontweight='bold') # Título Subgrafo con Zoom.
##------------------------------------------------------------------------------------------------------------------##
## 15º) ETIQUETAS de las ARISTAS:
    # a) Coordenadas de los nodos:
    node_coords= {node: (G.nodes[node]['x'], G.nodes[node]['y']) for node in Subgrafo_Zoom.nodes}

    # b) Etiquetas únicas (Evitar duplicidades):
    etiquetas_unicas_b= set()
    labels= []
    for u, v, clave, datos in Subgrafo_Zoom.edges(keys=True, data=True):
        if 'name' in datos:               # Si la arista tiene el atributo 'name'...
            nombre_arista= datos['name']  # Asígnalo como nombre de la arista.
            if str(nombre_arista) not in etiquetas_unicas_b: # Si NO está ya en la lista de aristas--> Añádelo.
                etiquetas_unicas_b.add(str(nombre_arista))
                x, y= (node_coords[u][0] + node_coords[v][0])/2, (node_coords[u][1] + node_coords[v][1])/2 # Posición de los nombres.
                labels.append(plt.text(x, y, nombre_arista, fontsize=8, color='black', ha='center', va='center'))

    # c) Ajustar automáticamente las etiquetas de las aristas para evitar superposiciones:
    adjust_text(labels, ax=ax, autoalign='xy', only_move={'text': 'xy'})
    st.pyplot(fig)  # Mostrar el ploteo en la APP WEB.

## GUARDAR IMAGEN:
    if Guardar_img==True:
            plt.savefig('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/Grafo_{}_Nodo_Mayor_Centralidad_ZOOM.png'.format(Nombre_Lugar, Nombre_Lugar))
##=========================================================================================================================================##
##=========================================================================================================================================##

### F) GRAFO CON LOS NODOS SEGÚN SU GRADO DE CENTRALIDAD MEDIA:
## 16º) DIBUJAR: Grafo con los Nodos según la Centralidad Media:
    fig, ax = ox.plot_graph(G,                                                # Grafo.
                            figsize=(Tamaño_figura[0]+15, Tamaño_figura[1]+15), # Tamaño figura.
                            node_color=ox.plot.get_node_colors_by_attr(G, 'Media_Centralidad', cmap='coolwarm'), # Color nodos.
                            node_size=Tamaño_nodo-20,                         # Tamaño nodos.
                            show=False, close=False, bgcolor='whitesmoke')    # Color fondo.
    plt.title('Nodos según su Centralidad Media en {}'.format(Nombre_Lugar), fontweight='bold') # Título.
##------------------------------------------------------------------------------------------------------------------##
## 17º) Añadir una BARRA de COLOR:
    sm=plt.cm.ScalarMappable(cmap='coolwarm',norm=plt.Normalize(vmin=0,vmax=100)); sm.set_array([])
    cbar= plt.colorbar(sm, shrink=0.5, ax=ax); cbar.ax.set_ylabel('Centralidad Media', rotation=270, labelpad=15)
    st.pyplot(fig)  # Mostrar el ploteo en la APP WEB.

## GUARDAR IMAGEN:
    if Guardar_img==True:
            plt.savefig('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/Grafo_{}_Nodos_según_Centralidad_MEDIA.png'.format(Nombre_Lugar, Nombre_Lugar))
##------------------------------------------------------------------------------------------------------------------##
## 18º) Crear MAPA INTERACTIVO con los NODOS según su Grado de CENTRALIDAD MEDIA:     (Mapa_Color,     Leyenda,     Estilo_Bordes_Nodos,                       Estilo_Mapa,
    mapa_nodos_segun_centralidad = geomapa_nodos.explore(column='Media_Centralidad', cmap='coolwarm', legend=True, style_kwds={'color':'black','weight':0.3}, tiles='CartoDB positron',
                                                     popup=['Nº Node', 'Media_Centralidad'], tooltip=None) # Información Nodos Mostrar, NO Mostrar resto de Info).
    
    # Guardar el Mapa Interactivo como archivo HTML (Google):
    mapa_nodos_segun_centralidad.save('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/MAPAS_INTERACTIVOS/MAPA_INTERACTIVO_Nodos_según_Centralidad_MEDIA.html'.format(Nombre_Lugar))
##=========================================================================================================================================##
##=========================================================================================================================================##

### G) GRAFO LLEGADA EN DISTINTOS MINUTOS DE TIEMPO:
## 19º) Añado el ATRIBUTO 'Tiempo_Andando' a las Aristas:
    meters_per_minute = Velocidad * 1000 / 60  # km per hour to m per minute
    for _, _, _, data in G.edges(data=True, keys=True):
        data["Tiempo_Andando"] = data["length"] / meters_per_minute
##------------------------------------------------------------------------------------------------------------------##
## 20º) FUNCIÓN PARA CREAR LOS POLÍGONOS DE LAS ISOCRONAS:

    ## 0º) Importar Librerías:
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path
    import geopandas as gpd
    from shapely.geometry import Point
    from shapely.geometry import LineString
    
    def Polígonos_Isocronas(G, Nodo_Clave, edge_buff=25, node_buff=50, infill=False):
        isochrone_polys = []  # Lista vacía donde irán los POLÍGONOS de las ISOCRONAS.
        for trip_time in sorted(Tiempos_Viaje, reverse=True):  # Para cada Tiempo de Viaje...

            ## 1º) SUBGRAFO VECINOS:          # ATRIBUTO (Aristas): 'Tiempo_Andando'.
            subgraph = nx.ego_graph(G, Nodo_Clave, radius=trip_time, distance="Tiempo_Andando")

            ## 2º) Obtención de las Coordenadas de los Nodos y Aristas del Subgrafo:
            node_points = [Point((data["x"], data["y"])) for node, data in subgraph.nodes(data=True)]
            nodes_gdf = gpd.GeoDataFrame({"id": list(subgraph.nodes)}, geometry=node_points)
            nodes_gdf = nodes_gdf.set_index("id")
            edge_lines = []
            for n_fr, n_to in subgraph.edges():
                f = nodes_gdf.loc[n_fr].geometry
                t = nodes_gdf.loc[n_to].geometry
                edge_lookup = G.get_edge_data(n_fr, n_to)[0].get("geometry", LineString([f, t]))
                edge_lines.append(edge_lookup)

            ## 3º) BUFFER (Área alrededor) de los Nodos y Aristas:
            n = nodes_gdf.buffer(node_buff).geometry
            e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
            all_gs = list(n) + list(e)
            new_iso = gpd.GeoSeries(all_gs).unary_union

            ## 4º) Rellenar los Polígonos (SIN Huecos):
            if infill:
                if new_iso is not None:
                    new_iso = [new_iso] if new_iso.geom_type == 'Polygon' else list(new_iso)
                else:
                    new_iso = []
            else:
                if new_iso is not None:
                    new_iso = [new_iso]
                else:
                    new_iso = []

            isochrone_polys.extend(new_iso)
        return isochrone_polys

##------------------------------------------------------------------------------------------------------------------##
## 21º) DIBUJO el GRAFO:
    fig, ax= ox.plot_graph(G,                  # Grafo.
                        figsize=Tamaño_figura, # Tamaño figura.
                        edge_color="silver",# Color Aristas.
                        edge_alpha=0.25,    # Transparencia Aristas.
                        node_color=['r' if node==Nodo_Clave else 'snow' for node in G.nodes], # Color nodos.
                        node_size=[8 if node==Nodo_Clave else 1 for node in G.nodes],         # Tamaño Nodos.
                        show=False, close=False, bgcolor='whitesmoke')                        # Color fondo.
    Tiempos= ", ".join(map(str, Tiempos_Viaje))                                               # String con los Tiempos de Viaje: a, b, c...
    plt.title('Puntos de llegada en {} en {} minutos'.format(Nombre_Lugar,Tiempos), fontweight='bold') # Título.
##------------------------------------------------------------------------------------------------------------------##
## 22º) Obtengo 1 COLOR para cada ISOCRONA:
    Iso_colores= ox.plot.get_colors(n=len(Tiempos_Viaje), cmap="plasma", start=0, return_hex=True)
##------------------------------------------------------------------------------------------------------------------##
## 23º) Obtengo los POLÍGONOS de las ISOCRONAS (con la función anterior):
    ISOCRONAS_Polig= Polígonos_Isocronas(G,                       # Grafo.
                                        Nodo_Clave= Nodo_Clave,# Nodo Salida.
                                        edge_buff=25, node_buff=0, infill=True) # Buffers.
##------------------------------------------------------------------------------------------------------------------##
## 24º) DIBUJO las ISOCRONAS:
    for polygon, fc in zip(ISOCRONAS_Polig, Iso_colores):
        if polygon.geom_type=='Polygon':
            coords= list(polygon.exterior.coords)
            path= Path(coords)                      # Transparencia.
            patch= PathPatch(path, fc=fc, ec="none", alpha=0.9, zorder=-1)
            ax.add_patch(patch)
        elif polygon.geom_type=='MultiPolygon':
            for p in polygon:
                coords= list(p.exterior.coords)
                path= Path(coords)                  # Transparencia.
                patch= PathPatch(path, fc=fc, ec="none", alpha=0.9, zorder=-1)
                ax.add_patch(patch)
##------------------------------------------------------------------------------------------------------------------##
## 25º) Añado una BARRA de COLOR:
    norm= plt.Normalize(min(Tiempos_Viaje), max(Tiempos_Viaje))
    sm=plt.cm.ScalarMappable(norm=norm, cmap='plasma_r'); sm.set_array([])
    cbar= plt.colorbar(sm, shrink=0.5, ax=ax); cbar.ax.set_ylabel('Tiempo de Viaje (minutos)', rotation=270, labelpad=15)
    st.pyplot(fig)  # Mostrar el ploteo en la APP WEB.

## GUARDAR IMAGEN:
    Tiempos= "_".join(map(str, Tiempos_Viaje))                                                # String con los Tiempos de Viaje: a_b_c...
    if Guardar_img==True:
            plt.savefig('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/Grafo_{}_Llegada_en_Tiempos_{}_minutos.png'.format(Nombre_Lugar, Nombre_Lugar, Tiempos))
##------------------------------------------------------------------------------------------------------------------##
## 26º) MAPA INTERACTIVO ISOCRONAS (Llegada en un determinado TIEMPO):
    # a) Obtener el color del nodo según la Isocrona correspondiente:
    node_colors = {}
    for trip_time, color in zip(sorted(Tiempos_Viaje, reverse=True), Iso_colores):
        subgraph = nx.ego_graph(G, Nodo_Clave, radius=trip_time, distance="Tiempo_Andando")
        for node in subgraph.nodes():
            node_colors[node] = color

    # b) Agregar nodos sin color definido con 'none':
    for node in geomapa_nodos.index:
        if node not in node_colors:
            node_colors[node] = 'none'
    node_colors[Nodo_Clave] = 'red' # Color del Nodo_Clave.

    # c) Asignar los colores a la columna--> 'Colores_Isocronas':
    geomapa_nodos['Colores_Isocronas']= node_colors

    # d) Función para asignar el TIEMPO en la lista Tiempos_Viaje a los colores:
    def asignar_indice_por_color(color):
        if color == 'none':
            return 'No se llega en estos tiempos'
        elif color == 'red':
            return 'Nodo de Salida'
        else:
            return Tiempos_Viaje[::-1][Iso_colores.index(color)]
    geomapa_nodos['Tiempo de Llegada (min)'] = geomapa_nodos['Colores_Isocronas'].apply(asignar_indice_por_color)

    # e) MAPA INTERACTIVO:           (Color_Nodos,               Leyenda,     Estilo_Borde_Nodos,        Estilo del Mapa,          Información al Pulsar,                  Resto de Info=NO):
    Mapa_ISOCRONAS= geomapa_nodos.explore(color='Colores_Isocronas', legend=True, style_kwds={'color':None}, tiles='CartoDB positron', popup=['Nº Node', 'Tiempo de Llegada (min)'], tooltip=None)

    # f) Guardar el Mapa Interactivo como archivo HTML (Google):
    Mapa_ISOCRONAS.save('./ANALISIS CENTRALIDAD/ANALISIS_CENTRALIDAD_{}/MAPAS_INTERACTIVOS/MAPA_INTERACTIVO_Llegada_en_Tiempos_{}_minutos.html'.format(Nombre_Lugar, Tiempos))
##===========================================================================================================================================##
##===========================================================================================================================================##

#### ESTILO DE TEXTOS Y BOTONES ####
### Color del FONDO del BOTÓN: ##
st.markdown("""<style>div.stButton > button:first-child {background-color: lightblue;}</style>""", unsafe_allow_html=True)

### Color y Centrado del TÍTULO a 'navy':                                  (La cosa está en el h1 [porque es el 1º texto]):
st.markdown("""<style> /* Cambiar el color del título a azul */ .element-container .stMarkdown h1 {color: navy; text-align: center;}</style>""",
    unsafe_allow_html=True)

### Color del 'header' a 'navy':                                             (La cosa está en el h2 [porque es el 2º texto]):
st.markdown("""<style>/* Cambiar el color del encabezado a azul */.element-container .stMarkdown h2 {color: navy;}</style>""", unsafe_allow_html=True)
##===========================================================================================================================================##
##===========================================================================================================================================##

def main():
    st.title('ANÁLISIS DE CENTRALIDAD Y DISTANCIA') # TÍTULO.
    st.image('https://raw.githubusercontent.com/Miguelgargor/apliacaion/main/imagen_mapa.png', use_column_width=True) # IMAGEN.

## DESCRIPCIÓN APP:
    st.write('Con esta aplicación, podrás seleccionar cualquier ciudad, pueblo o distrito en el Mundo y, al elegir los parámetros deseados, obtendrás:',
             '\n- El grafo descargado en caso de necesitarlo a posteriori.',
             '\n- Un mapa marcado con el punto geográfico desde el que se puede acceder a cualquier otro lugar de la zona de la manera más sencilla, reduciendo el tiempo y la distancia requeridos para llegar al destino, así como su ampliación.',
             '\n- Otro mapa con todos los puntos geográficos marcados según su grado de importancia en el estudio.',
             '\n- Y otro indicando los lugares alcanzables desde el punto marcado en diferentes tiempos y a una velocidad especifica.',
             '\n\n\nAdemás, se descargarán todos estos mapas en una carpeta en tu escritorio, que incluirá su versión interactiva.')

    st.header('Elige tus propios parámetros:') # SUBTÍTULO.

## Entrada de parámetros:
    # a) Lugar:
    lugar = st.text_input('**1- Lugar a estudiar**', value= 'Ávila, Castilla y León, España')
    # b) Tiempos de viaje: (Por defecto: [5,10,15]):
    tiempos_viaje = st.multiselect('**2- Tiempos de Viaje** (¿dónde llega desde el punto principal?)', list(range(121)), default=[5, 10, 15])
    tiempos_viaje = sorted(tiempos_viaje) # Ordenar los Tiempos de Viaje (ascendente).
    # c) Velocidad andando: (Por defecto: 6 Km/h):
    velocidad = st.number_input('**3- Velocidad de Viaje** (en Km/h)',min_value=1.0, max_value=100.0, value=6.0, step=0.1)
    # d) Distancia del cuadro del Grafo (compacto):
    distancia = st.slider('**4- Distancia máxima del lugar desde su punto medio** (en metros)', min_value=1000, max_value=8000, value=3500, step=100)
    # e) Tiempo máximo de espera al cargar el Grafo:
    tiempo_espera = st.slider('**5- Tiempo máximo de espera** (en segundos)', min_value=0, max_value=300, value=180, step=10)
    # f) Tipo de calles:
    tipo_calles = st.selectbox('**6- Tipo de calles**', ['drive', 'walk', 'bike', 'all', 'all_private'], index=0)
    # g) Tamaño del Nodo principal:    
    tamaño_nodo = st.slider('**7- Tamaño de los nodos en el grafo**', min_value=10, max_value=100, value=40, step=1)
    # h) Tamaño de la figura:
    tamaño_figura = st.selectbox('**8- Tamaño de la figura**', [(i, i) for i in range(1, 17)], index=9)
    # i) Guardar los Ploteos:
    guardar_img = True if st.radio('**9- ¿Guardar mapas?**', ['SI', 'NO'], index=0) == 'SI' else False
    # j) Guardar el Grafo:
    guardar_grafo = True if st.radio('**10- ¿Guardar Grafo?**', ['SI', 'NO'], index=0) =='SI' else False

### Crear el Botón para ejecutar el Código:
    if st.button('**GENERAR RESULTADOS**'):
        mensaje = st.empty()  # Crear un espacio en blanco para el mensaje que aparecerá.
        # Mensaje de CARGANDO...:
        mensaje.write('<p style="font-size: 20px; color: red; font-weight: bold;">Obtener los resultados puede tardar unos minutos...</p>',
                        unsafe_allow_html=True)
#### CARGADO DE LA FUNCIÓN EN SÍ:
        Dibujar_Nodo_Central(lugar, tiempos_viaje, velocidad, distancia, tiempo_espera, tipo_calles,
                             tamaño_nodo, tamaño_figura, guardar_img, guardar_grafo)
        
    ## Mensaje tras haber obtenido los resultados: (Indicando TAMAÑO, COLOR (limegreen), NEGRITA...):
        mensaje.write('<p style="font-size: 18px; color: limegreen; font-weight: bold;">¡Resultados obtenidos!</p>'
                    '<p style="font-size: 15px; color: limegreen;">(Puede ver los <strong>MAPAS INTERACTIVOS</strong> en la carpeta creada en tu escritorio)</p>',
                    unsafe_allow_html=True)
if __name__ == "__main__":
    main()
