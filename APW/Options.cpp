#include "Options.h"

#include <wx/stdpaths.h> 


Options::Options()
	: nrThreads(4), nrPoints(600), pathNo(10),
	paths{ { 
			{"K", "W", "X", "G", "L", "W"}, 
			{"W", "G", "X", "W", "L", "G"}, 
			{"W", "L", "G", "X", "W", "K"},
			{"L", "G", "X", "K", "G" },
			{"L", "G", "X", "W", "K", "G"},
			{"L", "G", "X", "U", "K", "G"}, 
			{"L", "G", "X", "U", "G"},
			{"L", "K", "W", "G", "X", "W", "L", "G", "K"},
			{"G", "X", "K", "G", "L", "K", "W", "X"},
			{"G", "X", "W", "L", "G", "K", "W", "U"},
			{"G", "X", "W", "L", "G", "K"},
			{"G", "X", "W", "L", "G", "K", "X"},
			{"G", "X", "W", "L", "G", "X"},
			{"G", "X", "W", "G", "U", "X"},
			{"G", "X", "W", "K", "L", "G"},
			{"G", "X", "W", "K", "G", "L", "U", "W", "L", "K"},
			{"G", "X", "U", "K", "G", "L", "W", "X"},
		} },
	method(1)
{
}

void Options::Open()
{
	if (m_fileconfig) return;

	const wxString dir = wxStandardPaths::Get().GetConfigDir() + wxFileName::GetPathSeparator();

	if(!wxFileName::DirExists(dir))
		wxFileName::Mkdir(dir, 0777, wxPATH_MKDIR_FULL);

	const wxString iniFilePath = dir + "APW.ini";

	m_fileconfig = new wxFileConfig("APW", wxEmptyString, iniFilePath);

	wxConfigBase::Set(m_fileconfig);
}

void Options::Close()
{
	delete m_fileconfig;
	m_fileconfig = nullptr;
	wxConfigBase::Set(nullptr);
}

void Options::Load()
{
	Open();
	wxConfigBase *conf=wxConfigBase::Get(false);
	if (conf)
	{
		nrThreads = conf->ReadLong("/nrThreads", 4);
		nrPoints = conf->ReadLong("/nrPoints", 600);
		pathNo = conf->ReadLong("/pathNo", 10);
		method = conf->ReadLong("/method", 1);
	}
	Close();
}

void Options::Save()
{
	Open();
	wxConfigBase *conf=wxConfigBase::Get(false);
	if (conf)
	{
		conf->Write("/nrThreads", static_cast<long int>(nrThreads));
		conf->Write("/nrPoints", static_cast<long int>(nrPoints));
		conf->Write("/pathNo", static_cast<long int>(pathNo));
		conf->Write("/method", static_cast<long int>(method));
	}

	if (m_fileconfig)
		m_fileconfig->Flush();
	Close();
}
