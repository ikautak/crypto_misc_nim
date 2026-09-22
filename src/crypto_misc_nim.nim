# This is just an example to get you started. A typical hybrid package
# uses this file as the main entry point of the application.

import crypto_misc_nim/submodule

when isMainModule:
  echo(getWelcomeMessage())
  var s = stdin.readLine()
  echo("You entered: ", s)
