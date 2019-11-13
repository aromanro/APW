#include "Options.h"

#include <wx/stdpaths.h> 


Options::Options()
	: nrThreads(4), nrPoints(600), pathNo(7),
	paths{ { {"K", "W", "X", "G", "L", "W"}, 
			 {"L", "G", "X", "K", "G" }, 
			 {"W", "G", "X", "W", "L", "G"}, 
			 {"L", "G", "X", "W", "K", "G"}, 
			 {"L", "G", "X", "U", "K", "G"}, 
			 {"G", "X", "K", "G", "L", "K", "W", "X"}, 
			 {"G", "X", "W", "L", "G", "K", "W", "U"},
			 {"G", "X", "W", "L", "G", "K"},
			 {"G", "X", "W", "G", "U", "X"}
		} },
	m_fileconfig(nullptr)
{
}

void Options::Open()
{
	if (m_fileconfig) return;

	wxString dir = wxStandardPaths::Get().GetConfigDir() + wxFileName::GetPathSeparator();

	if(!wxFileName::DirExists(dir))
		wxFileName::Mkdir(dir, 0777, wxPATH_MKDIR_FULL);

	wxString iniFilePath = dir + "APW.ini";

	m_fileconfig = new wxFileConfig("APW", wxEmptyString, iniFilePath);

	wxConfigBase::Set(m_fileconfig);
}

void Options::Close()
{
	delete m_fileconfig;
	m_fileconfig = NULL;
	wxConfigBase::Set(NULL);
}

void Options::Load()
{
	wxConfigBase *conf=wxConfigBase::Get(false);
	if (conf)
	{
		nrThreads = conf->ReadLong("/nrThreads", 20);
		nrPoints = conf->ReadLong("/nrPoints", 600);
		pathNo = conf->ReadLong("/pathNo", 7);
	}
}

void Options::Save()
{
	wxConfigBase *conf=wxConfigBase::Get(false);
	if (conf)
	{
		conf->Write("/nrThreads", static_cast<long int>(nrThreads));
		conf->Write("/nrPoints", static_cast<long int>(nrPoints));
		conf->Write("/pathNo", static_cast<long int>(pathNo));
	}

	if (m_fileconfig)
		m_fileconfig->Flush();
}
