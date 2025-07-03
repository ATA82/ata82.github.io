# _plugins/mathjax_injector.rb
module MathJaxInjector
  def self.inject_mathjax(content)
    # Check if there's a closing head tag, and insert the MathJax scripts before it.
    mathjax_code = <<~HTML
      <!-- MathJax Configuration -->
      <script type="text/x-mathjax-config">
        MathJax = {
          tex: {
            inlineMath: [['$', '$'], ['\\(', '\\)']]
          }
        };
      </script>
      <!-- MathJax Script -->
      <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    HTML

    # Insert the code before </head>, if present.
    if content.include?('</head>')
      content.gsub('</head>', "#{mathjax_code}\n</head>")
    else
      # Otherwise, simply append at the end.
      content + mathjax_code
    end
  end
end

# Use a hook to process pages after they've been rendered.
Jekyll::Hooks.register [:pages, :posts], :post_render do |page|
  if page.output_ext == '.html'
    page.output = MathJaxInjector.inject_mathjax(page.output)
  end
end

