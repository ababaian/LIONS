#! /bin/awk -f

BEGIN{
  FS=OFS="\t";
  t = 0
}; 

{ 
  s = $1; 
  if (s in a){
    t++
  } else {
    t=1
  }; 
  $1 = $1"_"t; 
  a[s][t] = $0; 
}; 

END{
  for (i in a) {
    for ( j in a[i]) {
      print a[i][j]
    }
  }
};
