# Package

version       = "0.1.0"
author        = "Katsuaki Oshio"
description   = "A new awesome nimble package"
license       = "MIT"
srcDir        = "src"
installExt    = @["nim"]
bin           = @["crypto_misc_nim", "miller_rabin_primality_test", "pollard_p_1", "lll"]


# Dependencies

requires "nim >= 2.2.12"
requires "bigints >= 1.1.0"
