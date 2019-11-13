#include "APWApp.h"
#include "APWFrame.h"


wxIMPLEMENT_APP(APWApp);


bool APWApp::OnInit()
{
	if (!wxApp::OnInit())
		return false;

	frame = new APWFrame("APW", wxPoint(50, 50), wxSize(1024, 800));
	frame->Show(true);

	return true;
}

