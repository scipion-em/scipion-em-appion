import os
from pyworkflow.plugin import Plugin
from bibtex import _bibtex


DOGPICKER_HOME_VAR = "DOGPICKER_HOME"

class AppionPlugin(Plugin):
    """ Appion plugin """

    def __init__(self):
        # Default variable definition
        configVars = {DOGPICKER_HOME_VAR: os.path.join(os.environ['EM_ROOT'],
                                                       'dogpicker-0.2.1')}

        Plugin.__init__(self, 'appion',
                              version="0.1",
                              logo="appion_logo.png",
                              configVars=configVars,
                              bibtex=_bibtex
                        )
    def registerPluginBinaries(self, env):

        # Add dogpicker
        env.addPackage('dogpicker', version='0.2.1',
                       tar='dogpicker-0.2.1.tgz',
                       default=True)


# Declare the plugin.
_plugin = AppionPlugin()
