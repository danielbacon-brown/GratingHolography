// -w320 -h240

#version 3.6;

#include "colors.inc"
#include "textures.inc"
#include "shapes.inc"

global_settings {max_trace_level 5 assumed_gamma 1.0}

camera {
	location <-1.600205, 3.200411, -4.800616>
	direction <0, 0,  2.25>
	right x*1.33
	look_at <0,0,0>
}

#declare Dist=80.0;
light_source {< -25, 50, -50> color White
	fade_distance Dist fade_power 2
}
light_source {< 50, 10,  -4> color Gray30
	fade_distance Dist fade_power 2
}
light_source {< 0, 100,  0> color Gray30
	fade_distance Dist fade_power 2
}

sky_sphere {
	pigment {
		gradient y
		color_map {
			[0, 1  color White color White]
		}
	}
}

#declare Xaxis = union{
	cylinder{
		<0,0,0>,<0.8,0,0>,0.05
	}
	cone{
		<0.8,0,0>, 0.1, <1,0,0>, 0
	}
	texture { pigment { color Red } }
}
#declare Yaxis = union{
	cylinder{
		<0,0,0>,<0,0.8,0>,0.05
	}
	cone{
		<0,0.8,0>, 0.1, <0,1,0>, 0
	}
	texture { pigment { color Green } }
}
#declare Zaxis = union{
	cylinder{
	<0,0,0>,<0,0,0.8>,0.05
	}
	cone{
		<0,0,0.8>, 0.1, <0,0,1>, 0
	}
	texture { pigment { color Blue } }
}
#declare Axes = union{
	object { Xaxis }
	object { Yaxis }
	object { Zaxis }
}
#declare Material_Vacuum = texture{ pigment{ color transmit 1.0 } }
#declare Material_Glass = texture{ pigment{ rgb <0.193304,0.563585,0.001251> } }
#declare Material_ITO = texture{ pigment{ rgb <0.479873,0.585009,0.808741> } }
#declare Material_SU8 = texture{ pigment{ rgb <0.822840,0.895962,0.350291> } }
#declare Material_PDMS = texture{ pigment{ rgb <0.858943,0.174108,0.746605> } }
#declare Material_Graphene = texture{ pigment{ rgb <0.303995,0.513535,0.710501> } }
#declare Layer_Front = union{
/*
	difference{
		intersection{
			plane{ <0.533400,0.000000,0>, 0.266700 }
			plane{ <-0.533400,-0.000000,0>, 0.266700 }
			plane{ <0.266700,0.461940,0>, 0.266701 }
			plane{ <-0.266700,-0.461940,0>, 0.266701 }
			plane{ <0.266700,-0.461940,0>, 0.266701 }
			plane{ <-0.266700,0.461940,0>, 0.266701 }
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.000000 }
		}
// nshapes = 0
		texture { Material_Vacuum }
	}
*/
	translate +z*0.000000
}
#declare Layer_Grating = union{
/*
	difference{
		intersection{
			plane{ <0.533400,0.000000,0>, 0.266700 }
			plane{ <-0.533400,-0.000000,0>, 0.266700 }
			plane{ <0.266700,0.461940,0>, 0.266701 }
			plane{ <-0.266700,-0.461940,0>, 0.266701 }
			plane{ <0.266700,-0.461940,0>, 0.266701 }
			plane{ <-0.266700,0.461940,0>, 0.266701 }
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.064000 }
		}
// nshapes = 1
prism{
	linear_spline
	0, 0.064000, 10,
	<0.000000,0.000000>,
	<0.533400,0.000000>,
	<0.533400,0.125509>,
	<0.533400,0.125509>,
	<0.236695,0.125529>,
	<0.236695,0.125529>,
	<0.236695,0.295100>,
	<0.236695,0.295100>,
	<0.000000,0.295379>,
	<0.000000,0.295379> 
	rotate +z*0.000000
	translate +x*-0.266700
	translate +y*-0.230900
}
		texture { Material_Vacuum }
	}
*/
	difference{
		intersection{
prism{
	linear_spline
	0, 0.064000, 10,
	<0.000000,0.000000>,
	<0.533400,0.000000>,
	<0.533400,0.125509>,
	<0.533400,0.125509>,
	<0.236695,0.125529>,
	<0.236695,0.125529>,
	<0.236695,0.295100>,
	<0.236695,0.295100>,
	<0.000000,0.295379>,
	<0.000000,0.295379> 
	rotate +z*0.000000
	translate +x*-0.266700
	translate +y*-0.230900
}
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.064000 }
		}
		texture { Material_SU8 }
	}
	translate +z*0.000000
}
#declare Layer_PrInterference = union{
	difference{
		intersection{
			plane{ <0.533400,0.000000,0>, 0.266700 }
			plane{ <-0.533400,-0.000000,0>, 0.266700 }
			plane{ <0.266700,0.461940,0>, 0.266701 }
			plane{ <-0.266700,-0.461940,0>, 0.266701 }
			plane{ <0.266700,-0.461940,0>, 0.266701 }
			plane{ <-0.266700,0.461940,0>, 0.266701 }
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 11.875000 }
		}
// nshapes = 0
		texture { Material_SU8 }
	}
	translate +z*0.064000
}
#declare Layer_Back = union{
	difference{
		intersection{
			plane{ <0.533400,0.000000,0>, 0.266700 }
			plane{ <-0.533400,-0.000000,0>, 0.266700 }
			plane{ <0.266700,0.461940,0>, 0.266701 }
			plane{ <-0.266700,-0.461940,0>, 0.266701 }
			plane{ <0.266700,-0.461940,0>, 0.266701 }
			plane{ <-0.266700,0.461940,0>, 0.266701 }
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.000000 }
		}
// nshapes = 0
		texture { Material_SU8 }
	}
	translate +z*11.939000
}
#declare Layers = union {
	//object{ Layer_Front }
	object{ Layer_Grating }
	object{ Layer_PrInterference }
	//object{ Layer_Back }
}

Axes
Layers
