function traceScanVolume(xAxis,yAxis,zAxis,raster)
    import zaber.motion.Units;  
    cont = 1;
    while cont ~= 0
        disp('Tracing Scan Volume...')
        % move to first corner
        xAxis.moveAbsolute(raster.xlims(1), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(1), Units.LENGTH_MILLIMETRES)
        zAxis.moveAbsolute(raster.zlims(1), Units.LENGTH_MILLIMETRES)
        % Scan top surface
        xAxis.moveAbsolute(raster.xlims(2), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(2), Units.LENGTH_MILLIMETRES)
        xAxis.moveAbsolute(raster.xlims(1), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(1), Units.LENGTH_MILLIMETRES)
        % Scan bottom surface
        zAxis.moveAbsolute(raster.zlims(2), Units.LENGTH_MILLIMETRES)
        xAxis.moveAbsolute(raster.xlims(2), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(2), Units.LENGTH_MILLIMETRES)
        xAxis.moveAbsolute(raster.xlims(1), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(1), Units.LENGTH_MILLIMETRES)

        cont = input('Repeat? [Yes = 1, No = enter]:');
    end
end