double setYmin(TGraphErrors * g, TGraphErrors * gA, TGraphErrors * gB){
  double ymin = 1.;
  for (int i = 0; i<g->GetN(); i++) {
    if(ymin > g->GetY()[i]-g->GetEY()[i]) ymin = g->GetY()[i]-g->GetEY()[i];
    if(ymin > gA->GetY()[i]-gA->GetEY()[i]) ymin = gA->GetY()[i]-gA->GetEY()[i];
    if(ymin > gB->GetY()[i]-gB->GetEY()[i]) ymin = gB->GetY()[i]-gB->GetEY()[i];
  }
  if(ymin>0) {
    return 0.;
  } else if( ymin > -0.01 ){
    return -0.01;
  } else if( ymin > -0.1 ){
    return -0.1;
  } else if( ymin > -0.2 ){
    return -0.2;
  } else if( ymin > -1.0 ){
    return -1.0;
  }
  return -10.;
}
double setYmax(TGraphErrors * g, TGraphErrors * gA, TGraphErrors * gB){
  double ymax = -1;
  for (int i = 0; i<g->GetN(); i++) {
    if(ymax < g->GetY()[i]+g->GetEY()[i])   ymax = g->GetY()[i]  + g->GetEY()[i];
    if(ymax < gA->GetY()[i]+gA->GetEY()[i]) ymax = gA->GetY()[i] + gA->GetEY()[i];
    if(ymax < gB->GetY()[i]+gB->GetEY()[i]) ymax = gB->GetY()[i] + gB->GetEY()[i];
  }
  if(ymax<0) {
    return 0.;
  } else if( ymax < 0.01 ){
    return 0.01;
  } else if( ymax < 0.1 ){
    return 0.1;
  } else if( ymax < 0.2 ){
    return 0.2;
  } else if( ymax < 0.4 ){
    return 0.4;
  } else if( ymax < 0.6 ){
    return 0.6;
  } else if( ymax < 1.0 ){
    return 1.0;
  }
  return 10.;
}

double plotmin(double ymin) {
    if(ymin>0) {
      return 0;
    } else if (ymin > -0.0004) {
      return -0.0004;
    } else if (ymin > -0.001) {
      return -0.001;
    } else if (ymin > -0.008) {
      return -0.01;
    } else if (ymin > -0.02) {
      return -0.02;
    } else if (ymin > -0.04) {
      return -0.05;
    } else if (ymin > -0.08) {
      return -0.1;
    } else if (ymin > -0.18) {
      return -0.2;
    } else if (ymin > -0.3) {
      return -0.4;
    } else if (ymin > -0.4) {
      return -0.5;
    } else {
      return -1;
    }
}

double plotmax(double ymax) {
    if(ymax<0) {
      return 0;
    } else if (ymax < 0.0005) {
      return 0.0004;
    } else if (ymax < 0.001) {
      return 0.001;
    } else if (ymax < 0.002) {
      return 0.002;
    } else if (ymax < 0.0095) {
      return 0.01;
    } else if (ymax < 0.02) {
      return 0.02;
    } else if (ymax < 0.04) {
      return 0.05;
    } else if (ymax < 0.08) {
      return 0.1;
    } else if (ymax < 0.18) {
      return 0.2;
    } else if (ymax < 0.3) {
      return 0.4;
    } else if (ymax < 0.4) {
      return 0.5;
    } else if (ymax < 0.8) {
      return 1.0;
    } else if (ymax < 1.3) {
      return 1.5;
    } else if (ymax < 2.5) {
      return 3.0;
    } else if (ymax < 4.0) {
      return 5.0;
    } else {
      return 10;
    }
}
