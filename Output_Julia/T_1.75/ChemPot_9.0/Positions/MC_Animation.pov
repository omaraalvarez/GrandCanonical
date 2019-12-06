global_settings {
	ambient_light rgb <0.2, 0.2, 0.2>	max_trace_level 15
}

background { color rgb <1, 1, 1> }

#default { finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.7} }

camera {
	perspective
	location <0, -15.5, 0>
	look_at <0, 0, 0>
}

light_source {
	<0, -50.0, 0>
	color rgb <0.3, 0.3, 0.3>
	fade_distance 100.0
	fade_power 0
	parallel
	point_at <0, 0, 0>
}

light_source {
	<0, 50.0, 0>
	color rgb <0.3, 0.3, 0.3>
	fade_distance 100.0
	fade_power 0
	parallel
	point_at <0, 0, 0>
}

light_source {
	<-50.0, 0, 0>
	color rgb <0.3, 0.3, 0.3>
	fade_distance 100.0
	fade_power 0
	parallel
	point_at <0, 0, 0>
}

light_source {
	<50.0, 0, 0>
	color rgb <0.3, 0.3, 0.3>
	fade_distance 100.0
	fade_power 0
	parallel
	point_at <0, 0, 0>
}

#macro Particle(rx, ry, rz)
	intersection {
		sphere {
			<rx, ry, rz>, 0.5
			pigment {rgbt <107/255, 173/255, 197/255, 0> }
		}
		box {
			<-L/2, -L/2, L/2>,	<L/2, L/2, -L/2>
			pigment {rgbt <107/255, 173/255, 197/255, 0> }
		}
	}
#end

#macro Walls(L, h)
	union {
		triangle {
			<L / 2, L / 2, h / 2>, <-L / 2, L / 2, h / 2>, <L / 2, -L / 2, h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
		triangle {
			<-L / 2, -L / 2, h / 2>, <-L / 2, L / 2, h / 2>, <L / 2, -L / 2, h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
	}

        union{
		triangle {
			<L / 2, L / 2, -h / 2>, <-L / 2, L / 2, -h / 2>, <L / 2, -L / 2, -h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
		triangle {
			<-L / 2, -L / 2, -h / 2>, <-L / 2, L / 2, -h / 2>, <L / 2, -L / 2, -h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
	}
#end

#macro Box(L, h)
	union{
		triangle {
			<-L / 2, -L / 2, h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
		triangle {
			<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
	}
        union {
		triangle {
			<L / 2, -L / 2, h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
		triangle {
			<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
	}

	
        union {
		triangle {
			<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <-L / 2, L / 2, -h /2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
triangle {
			<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <L / 2, L / 2, h / 2>
			pigment { rgbt <0.75, 0.75, 0.75, 0.5> }
		}
	}
#end
#declare L = 10.0;
#declare h = 10.0;
#fopen File_Positions concat("/home/omaralvarez/git/Mezei/Output_Julia/T_1.75/ChemPot_9.0/Positions/Pos_", str(clock, 1, 0), ".xyz") read
	#while (defined( File_Positions ))
		#read (File_Positions, rx, ry, rz)
		Particle(rx, ry, -rz)
		#declare PBC = false;
		#if (rx > (L - 1) / 2)
			#declare rx = rx - L;
			#declare PBC = true;
		#end
		#if (rx < -(L - 1) / 2)
			#declare rx = rx + L;
			#declare PBC = true;
		#end
		
        #if (ry > (L - 1) / 2)
			#declare ry = ry - L;
			#declare PBC = true;
		#end
		#if (ry < -(L - 1) / 2)
			#declare ry = ry + L;
			#declare PBC = true;
		#end
		#if (rz > (h - 1) / 2)
			#declare rz = rx - h;
			#declare PBC = true;
		#end
		#if (rz < -(h - 1) / 2)
			#declare rz = rz + h;
			#declare PBC = true;
		#end
		#if (PBC)
			Particle(rx, ry, -rz)
		#end
	#end
#fclose File
Walls(L, h)
Box(L, h)
