; ==============================================================================
; Script: calc_WPSH_NBI.ncl
; 
; Description: This script calculates the Western Pacific Subtropical High 
;              Northern Boundary Index (WPSHNBI) based on ERA5 500-hPa 
;              Geopotential.
; 
; Input Data: Geopotential.nc (Variable: z)
; 
; Output:     WPSH_North_Boundary_Index.csv
; 
; Requirements: NCAR Command Language (NCL) Version 6.6.2 or later.
; ==============================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
    ; ---Read Data ---
    f_z = addfile("Geopotential.nc", "r")
    z = tofloat(f_z->z)
    
    lat = f_z->latitude
    lon = f_z->longitude
    
    z_ndims = dimsizes(dimsizes(z))
    z_gpm = z / tofloat(9.80665)
    
    ; ---Define Regions and Thresholds ---
    lat_min = 10.0
    lat_max = 60.0
    lon_min = 110.0
    lon_max = 150.0
    z_threshold = 5880.0
    
    lat_idx_int = tointeger(ind(lat.ge.lat_min .and. lat.le.lat_max))
    lon_idx_int = tointeger(ind(lon.ge.lon_min .and. lon.le.lon_max))
    
    ; ---Extract Regional Data ---
    if (z_ndims .eq. 3) then
        z_region = z_gpm(:, lat_idx_int, lon_idx_int)
    else if (z_ndims .eq. 4) then
        z_region = z_gpm(:, 0, lat_idx_int, lon_idx_int)
    end if
    end if
    
    lat_region = lat(lat_idx_int)
    lon_region = lon(lon_idx_int)
    
    ntime = tointeger(dimsizes(z_region(:,0,0)))
    nlat  = tointeger(dimsizes(lat_region))
    nlon  = tointeger(dimsizes(lon_region))
    
    ; --- Initialize Arrays ---
    north_boundary_lat = new((/ntime, nlon/), "float")
    north_boundary_lat@_FillValue = -999.0
    
    valid_north_count = new(ntime, "integer")
    has_5880_flag     = new(ntime, "integer")
    
    avg_north_boundary_lat = new(ntime, "float")
    avg_north_boundary_lat@_FillValue = -999.0
    
    print("Processing " + ntime + " time steps for WPSHNBI...")
    
    ; ---Main Calculation Loop ---
    do t = 0, ntime-1
        if (mod(t, 100) .eq. 0) then
            print("   Processing time step: " + t + " / " + ntime)
        end if
        
        valid_north_count(t) = 0
        has_5880_flag(t)     = 0
        sum_north_lat = 0.0
        
        ; Step 1: Find WPSH Northern Boundary
        do i = 0, nlon-1
            z_profile = z_region(t, :, i)
            n_crossings = 0
            crossing_lats = new(nlat-1, "float")
            crossing_lats@_FillValue = -999.0
            
            do j = 0, nlat-2
                ; 过滤缺测值
                if (ismissing(z_profile(j)) .or. ismissing(z_profile(j+1))) then
                    continue
                end if
                
                ; 过滤未跨越 5880 阈值的点
                if ((z_profile(j) - z_threshold) * (z_profile(j+1) - z_threshold) .gt. 0) then
                    continue
                end if
                
                ; 核心计算
                z1 = z_profile(j) - z_threshold
                z2 = z_threshold - z_profile(j+1)
                
                if (abs(z1) + abs(z2) .gt. 0) then
                    weight1 = z2 / (z1 + z2)
                    weight2 = z1 / (z1 + z2)
                    interp_lat = tofloat(lat_region(j) * weight1 + lat_region(j+1) * weight2)
                    
                    crossing_lats(n_crossings) = interp_lat
                    n_crossings = n_crossings + 1
                end if
            end do
            
            ; 保存最北的交点
            if (n_crossings .gt. 0) then
                has_5880_flag(t) = 1
                max_lat = max(crossing_lats(0:n_crossings-1))
                north_boundary_lat(t, i) = max_lat
                sum_north_lat = sum_north_lat + max_lat
                valid_north_count(t) = valid_north_count(t) + 1
            else
                north_boundary_lat(t, i) = north_boundary_lat@_FillValue
            end if
            delete(crossing_lats)
        end do
        
        ; Step 2: Calculate Averages
        if (has_5880_flag(t) .eq. 1 .and. valid_north_count(t) .gt. 0) then
            avg_north_boundary_lat(t) = sum_north_lat / tofloat(valid_north_count(t))
        else
            avg_north_boundary_lat(t) = avg_north_boundary_lat@_FillValue
        end if
    end do
    
    ; ---Time Formatting ---
    time_str = new(ntime, "string")
    if (isfilevar(f_z, "time")) then
        time_var = f_z->time
        if (isatt(time_var, "units")) then
            time_str = cd_string(time_var, "%Y-%m-%d")
        else
            do t = 0, ntime-1
                time_str(t) = "Step_" + tostring(t+1)
            end do
        end if
    else
        do t = 0, ntime-1
            time_str(t) = "Step_" + tostring(t+1)
        end do
    end if
    
    ; ---Output CSV Files ---
    north_filename = "WPSH_North_Boundary_Index.csv"
    if (fileexists(north_filename)) then
        system("rm -f " + north_filename)
    end if
    
    north_header = "Time,Avg_North_Boundary_Lat(N),Valid_North_Count,Has_5880"
    do i = 0, nlon-1
        north_header = north_header + ",Lon_" + sprintf("%.1f", lon_region(i)) + "_North_Boundary"
    end do
    write_table(north_filename, "w", [/north_header/], "%s")
    
    do t = 0, ntime-1
        avg_north_str = "NA"
        if (.not.ismissing(avg_north_boundary_lat(t))) then
            avg_north_str = sprintf("%.2f", avg_north_boundary_lat(t))
        end if
        north_line = time_str(t) + "," + avg_north_str + "," + tostring(valid_north_count(t)) + "," + tostring(has_5880_flag(t))
        do i = 0, nlon-1
            if (.not.ismissing(north_boundary_lat(t, i))) then
                north_line = north_line + "," + sprintf("%.2f", north_boundary_lat(t, i))
            else
                north_line = north_line + ",NA"
            end if
        end do
        write_table(north_filename, "a", [/north_line/], "%s")
    end do
    
    print("Processing complete. WPSH_North_Boundary_Index.csv generated.")
end
