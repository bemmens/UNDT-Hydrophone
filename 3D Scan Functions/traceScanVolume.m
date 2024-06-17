function traceScanVolume(xAxis,yAxis,zAxis,raster)
    import zaber.motion.Units;  
    cont = 1;
    while cont ~= 0
        disp('Tracing Scan Volume...')
        xAxis.moveAbsolute(raster.xlims(1), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(1), Units.LENGTH_MILLIMETRES)
        zAxis.moveAbsolute(raster.zlims(1), Units.LENGTH_MILLIMETRES)
        xAxis.moveAbsolute(raster.xlims(2), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(2), Units.LENGTH_MILLIMETRES)
        xAxis.moveAbsolute(raster.xlims(1), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(1), Units.LENGTH_MILLIMETRES)
    
        % Bottom done now top
        zAxis.moveAbsolute(raster.zlims(2), Units.LENGTH_MILLIMETRES)
        xAxis.moveAbsolute(raster.xlims(2), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(2), Units.LENGTH_MILLIMETRES)
        xAxis.moveAbsolute(raster.xlims(1), Units.LENGTH_MILLIMETRES)
        yAxis.moveAbsolute(raster.ylims(1), Units.LENGTH_MILLIMETRES)

        cont = input('Repeat? [Yes = 1, No = 0]:');
    end
end