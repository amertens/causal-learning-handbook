-- convert-interactive.lua
-- Converts R code blocks containing "#| interactive: true" into {webr-r}
-- blocks so the quarto-webr extension picks them up.

function CodeBlock(el)
  -- Only act on code blocks with r class or cell-code class
  local is_r = el.classes:includes("r") or el.classes:includes("cell-code")
  if not is_r then
    return nil
  end

  -- Check for #| interactive: true in the code text itself
  -- (Quarto processes #| options before Lua, but the text still contains them)
  if el.text:match("#|%s*interactive:%s*true") then
    -- Replace the class
    local new_classes = pandoc.List({})
    for _, c in ipairs(el.classes) do
      if c == "r" then
        new_classes:insert("webr-r")
      else
        new_classes:insert(c)
      end
    end
    el.classes = new_classes
    -- Remove the interactive line from code text
    el.text = el.text:gsub("#|%s*interactive:%s*true%s*\n", "")
    return el
  end

  return nil
end
