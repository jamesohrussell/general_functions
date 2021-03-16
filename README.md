# Generic Python functions for geophysical and data science

geophys_functions:
 * calc_area
 * calc_area_and_volrainrate
 * calc_distance
 * calc_direction
 * calc_distdir
 * calc_propagation
 * load_land
 * is_land
 * interp_TC
 * calc_if_TC
 * lonFlip

time_functions:
* time_since (replicates NCL cd_calendar)
* time_since_inv (replicates NCL cd_calendar_inv)
* utc_to_local - converts utc to local time using timezone
* calc_local_solar_time
* calc_local_solar_hour

shape_functions:
* label_wdiags
* find_corners
* points_in_shape
* fit_ellipse_svd_earth
* periodic_cmass_earth
* convert_cendirlen_latlon_earth

plot_functions
* truncate_colormap
* plot_pixel_axes_earth
* plot_pixel

misc_functions
* cartesian_distance
* cartesian_speed
* cartesian_direction
* divzero
* k_closest
* k_closest_ma
* create_2d_dataframe
* write_var
* write_var_compress
* write_group

